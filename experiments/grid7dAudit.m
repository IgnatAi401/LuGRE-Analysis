function report = grid7dAudit(windows, outputFolder)
% GRID7DAUDIT 批处理“位置+速度+钟差 7 维网格”DPE 的可行性与运算量审计（不改求解器）。
% 沿用 mainBatch 的观测、时间轴与动力学模型；以现有批处理 DPE（IRLS）解作参照点，回答：
%   1. 数据规模：历元、观测、每历元卫星数、样条钟差参数个数；
%   2. 钟差能否用 1 个参数表示：不同时长子窗口内钟差偏离常数/线性的幅度；
%   3. 网格点能否不重积分轨道：状态转移矩阵（STM）线性近似的误差；
%   4. 不同钟差模型下 6 维轨道状态的形式精度（可观测性）；
%   5. 初值偏差与所需网格范围，穷举与逐级细化的网格点数；
%   6. 单个网格点代价的实测耗时（每历元钟差用中位数精确消去）；
%   7. “6 维网格 + 每历元钟差”的 L1 最优解（网格收敛的终点）与真值的误差预览。
%      预览用 IRLS 求该凸问题的最优点，只作审计参照，不是拟议方法本身。
% 例：grid7dAudit({'OP5_0','OP9_0','OP14_0'})
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, genpath(fullfile(root,'src')));
if nargin < 1, windows = {'OP5_0','OP9_0','OP14_0'}; end
if nargin < 2, outputFolder = fullfile(root,'results','grid7d_audit'); end
if ~isfolder(outputFolder), mkdir(outputFolder); end
report = struct([]);
for w = 1:numel(windows)
    cfg = lugreConfig(); cfg.data.window = windows{w};
    cfg.run.plot = false; cfg.batch.methods = {'dpeIrls'};
    fprintf('\n===== %s =====\n', windows{w});
    timer = tic; base = runBatchPositioning(cfg); r.baselineSeconds = toc(timer);
    r.window = windows{w};
    m = numel(base.gpsSeconds); rotations = zeros(3,3,m);
    for k = 1:m, rotations(:,:,k) = icrfToItrfRotation(base.gpsSeconds(k), cfg.frame); end
    timer = tic;
    solveBatchIrls(base.observations, base.gpsSeconds, rotations, base.init.state, base.moon0, base.tRef, 'dpeIrls', cfg);
    r.irlsSolveSeconds = toc(timer);
    fprintf('[耗时] 现有 IRLS 批处理求解 %.1f s（整个 runBatchPositioning %.0f s）\n', r.irlsSolveSeconds, r.baselineSeconds);
    r.dataStats = dataScale(base);
    r.clock = clockStats(base);
    [r.stm, Phi, sol] = stmStats(base, cfg);
    [r.observability, H, y0, groups] = observability(base, rotations, cfg, Phi, sol);
    r.gridSize = gridSize(base, r.observability);
    r.timing = timing(H, y0, groups);
    r.perEpochOptimum = perEpochOptimum(base, rotations, cfg, Phi, H, y0, groups);
    if isfield(base,'truth'), r.truthMedianM = base.truth.dpeIrls.medianM; end
    report = [report, r]; %#ok<AGROW>
end
save(fullfile(outputFolder,'grid7d_audit.mat'),'report');
fprintf('\n结果已保存到 %s\n', fullfile(outputFolder,'grid7d_audit.mat'));
end

function s = dataScale(base)
counts = base.observationCount(:);
epochs = base.gpsSeconds(:);
observed = epochs(counts>0);
s.epochs = numel(epochs);
s.observations = sum(counts);
s.spanS = observed(end) - observed(1);
s.epochsBySatCount = accumarray(counts+1, 1)'; % 第 k 个元素：k-1 颗星的历元数
s.splineClockParameters = numel(base.dpeIrls.clockCoefficientsM);
fprintf('[规模] 时长 %.0f s，历元 %d，观测 %d，样条钟差参数 %d\n', s.spanS, s.epochs, ...
    s.observations, s.splineClockParameters);
fprintf('[规模] 每历元卫星数 0..%d 的历元数：%s\n', numel(s.epochsBySatCount)-1, mat2str(s.epochsBySatCount));
end

function s = clockStats(base)
% 以批处理 DPE 的样条钟差为“真实”钟差形状，看 1 个钟差参数能覆盖多长的子窗口。
t = base.gpsSeconds(:); c = base.dpeIrls.clockM(:);
good = isfinite(c); t = t(good); c = c(good);
lengths = [10 30 60 120 300 600 Inf];
s.lengthsS = lengths;
s.constantMaxM = nan(size(lengths)); s.linearMaxM = nan(size(lengths));
for j = 1:numel(lengths)
    L = min(lengths(j), t(end)-t(1)+1);
    edges = t(1):L:t(end)+L;
    maxC = []; maxL = [];
    for k = 1:numel(edges)-1
        in = t >= edges(k) & t < edges(k+1);
        if nnz(in) < 3, continue; end
        cc = c(in); tt = t(in) - mean(t(in));
        maxC(end+1) = max(abs(cc - median(cc))); %#ok<AGROW>
        p = [ones(size(tt)) tt] \ cc;
        maxL(end+1) = max(abs(cc - [ones(size(tt)) tt]*p)); %#ok<AGROW>
    end
    s.constantMaxM(j) = median(maxC); s.linearMaxM(j) = median(maxL);
end
s.fullRangeM = max(c) - min(c);
fprintf('[钟差] 全窗口变化 %.0f m\n', s.fullRangeM);
fprintf('[钟差] 子窗口长度/s        %s\n', sprintf('%8g', lengths));
fprintf('[钟差] 常数模型最大偏离/m  %s\n', sprintf('%8.1f', s.constantMaxM));
fprintf('[钟差] 线性模型最大偏离/m  %s\n', sprintf('%8.1f', s.linearMaxM));
end

function [s, Phi, sol] = stmStats(base, cfg)
% 有限差分 STM（参照解处），检验网格点用 r(t)=r_nom(t)+Phi(t)*dx 代替重积分的误差。
sol = base.dpeIrls; epochs = base.gpsSeconds(:); tRef = base.tRef; moon0 = base.moon0;
x0 = sol.state;
timer = tic; nominal = propagateArc(x0, moon0, tRef, epochs, cfg); s.propagateSeconds = toc(timer);
steps = [1e3*ones(3,1); ones(3,1)];
m = numel(epochs); Phi = zeros(3,6,m);
for j = 1:6
    dx = zeros(6,1); dx(j) = steps(j);
    p = propagateArc(x0+dx, moon0, tRef, epochs, cfg);
    Phi(:,j,:) = permute((p - nominal)/steps(j), [2 3 1]);
end
tau = epochs - tRef;
simple = zeros(3,6,m);
for k = 1:m, simple(:,:,k) = [eye(3), tau(k)*eye(3)]; end
radial = x0(1:3)/norm(x0(1:3));
tests = {[50e3*radial; 0;0;0], [0;0;0; 50*radial], [50e3*radial; 50*radial], [5e3*null(radial')*[1;1]/sqrt(2); 5*radial]};
labels = {'径向 50 km','径向速度 50 m/s','径向 50 km + 50 m/s','横向 5 km + 5 m/s'};
s.tests = labels; s.stmErrorM = nan(1,numel(tests)); s.simpleErrorM = s.stmErrorM;
for j = 1:numel(tests)
    truthArc = propagateArc(x0+tests{j}, moon0, tRef, epochs, cfg);
    predStm = nominal + squeeze(pagemtimes(Phi, tests{j}))';
    predSimple = nominal + squeeze(pagemtimes(simple, tests{j}))';
    s.stmErrorM(j) = max(vecnorm(truthArc - predStm,2,2));
    s.simpleErrorM(j) = max(vecnorm(truthArc - predSimple,2,2));
    fprintf('[STM] %-18s 最大位置误差：STM %.3f m，[I tau*I] %.1f m\n', labels{j}, ...
        s.stmErrorM(j), s.simpleErrorM(j));
end
fprintf('[STM] 单次整弧轨道积分 %.2f s（每个网格点重积分不可行）\n', s.propagateSeconds);
end

function [s, H, y0, groups] = observability(base, rotations, cfg, Phi, sol)
% 在参照解处线性化：伪距对 6 维初始状态的偏导 H；比较不同钟差模型下的形式精度。
obs = base.observations; used = sol.used;
epochs = base.gpsSeconds(:); tRef = base.tRef;
m = numel(epochs);
positionIcrf = propagateArc(sol.state, base.moon0, tRef, epochs, cfg);
idx = obs.epochIndex(used); sat = obs.satPosition(used,:); range = obs.range(used);
N = numel(range);
posEcef = zeros(m,3);
for k = 1:m, posEcef(k,:) = (rotations(:,:,k)*positionIcrf(k,:)')'; end
los = posEcef(idx,:) - sat; dist = vecnorm(los,2,2); unit = los./dist;
H = zeros(N,6);
for k = unique(idx)'
    rows = idx==k;
    H(rows,:) = unit(rows,:)*rotations(:,:,k)*Phi(:,:,k);
end
y0 = range - dist; % 参照解处残差（含钟差）
sigma = sol.robustSigmaM;
tau = epochs(idx) - tRef; halfSpan = max(abs(tau));
[~,~,groups] = unique(idx);
% 每历元自由钟差：减去历元内均值即消去钟差（仅 >=2 颗星的历元有信息）
Hc = H - groupMean(H, groups);
models = {'spline120s', H, splineBasis(tau, 120); ...
          'perEpoch', Hc, zeros(N,0); ...
          'constant', H, ones(N,1); ...
          'linear', H, [ones(N,1) tau/halfSpan]; ...
          'quadratic', H, [ones(N,1) tau/halfSpan (tau/halfSpan).^2]};
fprintf('[可观测] 稳健残差尺度 sigma=%.1f m；tRef 处 6 维状态形式标准差：\n', sigma);
fprintf('          %-11s %12s %12s %12s %12s\n','钟差模型','位置3D/m','最弱位置轴/m','最强位置轴/m','速度3D/(m/s)');
for j = 1:size(models,1)
    A = [models{j,2}, models{j,3}];
    n = vecnorm(A,2,1); n(n==0) = 1;
    C = sigma^2 * inv((A./n)'*(A./n)) ./ (n'*n);
    C = C(1:6,1:6);
    e = sqrt(sort(eig((C(1:3,1:3)+C(1:3,1:3)')/2),'descend'));
    s.(models{j,1}) = struct('covariance',C,'position3dM',sqrt(trace(C(1:3,1:3))), ...
        'positionAxesM',e','velocity3dMps',sqrt(trace(C(4:6,4:6))));
    fprintf('          %-11s %12.0f %12.0f %12.1f %12.3f\n', models{j,1}, s.(models{j,1}).position3dM, ...
        e(1), e(3), s.(models{j,1}).velocity3dMps);
end
s.sigmaM = sigma; s.N = N; s.Hc = Hc;
end

function g = groupMean(X, groups)
counts = accumarray(groups,1);
g = zeros(size(X));
for j = 1:size(X,2)
    mu = accumarray(groups, X(:,j))./counts;
    g(:,j) = mu(groups);
end
end

function B = splineBasis(tau, spacing)
span = [min(tau) max(tau)];
count = ceil((span(2)-span(1))/spacing) + 3;
u = (tau - span(1))/spacing + 3 - (0:count-1);
B = zeros(size(u)); piece = floor(u);
a = piece==0; B(a) = u(a).^3/6;
b = piece==1; B(b) = (-3*u(b).^3 + 12*u(b).^2 - 12*u(b) + 4)/6;
c = piece==2; B(c) = (3*u(c).^3 - 24*u(c).^2 + 60*u(c) - 44)/6;
d = piece==3; B(d) = (4-u(d)).^3/6;
B = B(:, any(B,1));
end

function s = gridSize(base, ob)
% 网格沿“每历元钟差”模型的法方程特征方向（白化坐标）布置。
% 范围：max(3 倍初值偏差, 5 sigma_i)；目标步长 0.5 sigma_i。
C = ob.perEpoch.covariance;                         % 6x6 形式协方差（m, m/s 混合单位）
[V, D] = eig((C+C')/2);
axisSigma = sqrt(diag(D));                          % 各主轴的形式标准差
dx0 = base.init.state - base.dpeIrls.state;
offsetSigma = abs(V'*dx0)./axisSigma;               % 初值偏差（以 sigma_i 为单位）
half = max(3*offsetSigma, 5);
exhaustive = prod(2*half/0.5 + 1);
levels = ceil(log2(max(half)/0.5/2)) + 1;
s.axisSigma = axisSigma'; s.initOffsetSigma = offsetSigma'; s.halfWidthSigma = half';
s.exhaustivePoints = exhaustive;
s.coarseToFine = struct('pointsPerAxis',5,'levels',levels,'points',levels*5^6);
% 物理坐标轴对齐网格（不白化）：位置 ±3|初值偏差|、步长 = 最小位置轴 sigma/2，速度同理
posHalf = max(3*abs(dx0(1:3)), 1e3); velHalf = max(3*abs(dx0(4:6)), 1);
posStep = 0.5*sqrt(min(eig(C(1:3,1:3)))); velStep = 0.5*sqrt(min(eig(C(4:6,4:6))));
s.axisAlignedPoints = prod(2*posHalf/posStep+1) * prod(2*velHalf/velStep+1);
fprintf('[网格] 主轴 sigma：        %s\n', sprintf('%10.3g', axisSigma));
fprintf('[网格] 初值偏差（sigma 单位）：%s\n', sprintf('%10.1f', offsetSigma));
fprintf('[网格] 白化后穷举 %.2e 点；5 点/轴逐级细化 %d 级 x 15625 = %.2e 点；不白化轴对齐穷举 %.2e 点\n', ...
    exhaustive, levels, s.coarseToFine.points, s.axisAlignedPoints);
end

function s = timing(H, y0, groups)
% 单点代价：r = y0 - H*dx；每历元钟差取残差中位数（L1 最优，等价于无限细的钟差网格）。
counts = accumarray(groups,1);
K = max(counts); mEff = numel(counts);
order = zeros(mEff, K);
[~, sorter] = sort(groups);
pos = 1;
for k = 1:mEff
    order(k,1:counts(k)) = sorter(pos:pos+counts(k)-1)'; pos = pos+counts(k);
end
pad = order == 0; order(pad) = numel(y0) + 1;
G = 2048;
dx = randn(6, G) .* [1e3*ones(3,1); ones(3,1)];
timer = tic;
R = [y0 - H*dx; nan(1,G)];
Rp = reshape(R(order,:), mEff, K, G);
med = median(Rp, 2, 'omitnan');
cost = squeeze(sum(abs(Rp - med), [1 2], 'omitnan')); %#ok<NASGU>
s.perPointUs = toc(timer)/G*1e6;
% 对照：线性化后固定钟差（真 7 维，常数钟差）只需 N 次运算
timer = tic;
R = y0 - H*dx - 5*randn(1,G);
cost = sum(abs(R),1); %#ok<NASGU>
s.perPoint7dUs = toc(timer)/G*1e6;
s.observations = numel(y0); s.epochGroups = mEff; s.maxSatsPerEpoch = K;
fprintf('[耗时] 每点：6 维+每历元中位数钟差 %.1f us；真 7 维（常数钟差）%.1f us（N=%d）\n', ...
    s.perPointUs, s.perPoint7dUs, numel(y0));
end

function s = perEpochOptimum(base, rotations, cfg, Phi, H, y0, groups)
% 线性化模型 y0 = H*dx + c_k + e 的 L1 最优（c_k 每历元自由）。问题是凸的，
% 网格逐级细化与 IRLS 收敛到同一最优点，这里用 IRLS 快速得到该点以预览精度。
dx = zeros(6,1); omega = ones(size(y0));
for it = 1:200
    wsum = accumarray(groups, omega);
    Hw = H - wmean(H); yw = y0 - wmean(y0);
    next = (sqrt(omega).*Hw) \ (sqrt(omega).*yw);
    r = yw - Hw*next;
    omega = 1./max(abs(r), 0.01);
    if norm(next - dx) < 1e-4, dx = next; break; end
    dx = next;
end
    function g = wmean(X)
        g = zeros(size(X));
        for j = 1:size(X,2)
            mu = accumarray(groups, omega.*X(:,j))./wsum; g(:,j) = mu(groups);
        end
    end
residual = y0 - H*dx; residual = residual - groupMedian(residual, groups);
s.cost = sum(abs(residual));
s.dxFromSpline = dx';
epochs = base.gpsSeconds(:); m = numel(epochs);
nominal = propagateArc(base.dpeIrls.state, base.moon0, base.tRef, epochs, cfg);
pos = zeros(m,3);
for k = 1:m, pos(k,:) = (rotations(:,:,k)*(nominal(k,:)' + Phi(:,:,k)*dx))'; end
s.positionEcef = pos;
fprintf('[预览] 每历元钟差 L1 最优相对样条解：位置 %.0f m，速度 %.3f m/s\n', norm(dx(1:3)), norm(dx(4:6)));
if isfield(base,'truth')
    t = base.truth.position; e = pos - t; good = all(isfinite(e),2);
    u = t./vecnorm(t,2,2); radial = sum(e.*u,2);
    s.truthMedianM = median(vecnorm(e(good,:),2,2));
    s.truthRadialRmsM = rms(radial(good));
    s.truthTransverseRmsM = rms(vecnorm(e(good,:) - radial(good).*u(good,:),2,2));
    fprintf('[预览] 真值误差中位数：每历元钟差 %.0f m（样条 IRLS %.0f m）；径向 RMS %.0f m，横向 RMS %.0f m\n', ...
        s.truthMedianM, base.truth.dpeIrls.medianM, s.truthRadialRmsM, s.truthTransverseRmsM);
end
end

function g = groupMedian(x, groups)
mu = accumarray(groups, x, [], @median); g = mu(groups);
end
