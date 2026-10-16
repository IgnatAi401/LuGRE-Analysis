function solution = solveBatchIrls(obs, epochs, rotations, state0, moon0, tRef, method, cfg, lunar)
% SOLVEBATCHIRLS 多历元批处理定轨：整段观测联合估计一组状态。
% 未知数：tRef 的地心 ICRF [r0;v0]、钟差 c*dt(t)（整段多项式或三次 B 样条）、
% 各信号相对参考信号的常值偏差。
% cfg.batch.dynamics = 'lunarStatic' 时状态改为月固系固定位置 p（3 维，m），位置
% x_k = lunar.moonEcef(k,:)' + lunar.moonFixedToEcef(:,:,k)*p（computeMoonFixedFrame），不积分；
% 此时 moon0 不使用，设计矩阵为稀疏（长窗口样条钟差有数百列），无观测支撑的样条列不参与。
% method = 'dpeIrls'：加权 L1 伪距残差代价（IRLS 求解），可选沿法方程最弱方向一维搜索（weakSearch）；
% method = 'ls' ：同一模型、同一观测的加权 L2 代价（高斯-牛顿），作为公平对照。
% obs 字段：epochIndex、satPosition(ECEF, 已含 Sagnac)、range(已加卫星钟差)、weight、signalId。
% rotations(:,:,k)：epochs(k) 的 ICRF->ECEF 矩阵。单位 m、m/s、s。
opts = cfg.batch;
assert(any(strcmp(method,{'dpeIrls','ls'})),'LuGRE:BatchMethod','method 必须为 dpeIrls 或 ls。');
static = isfield(opts,'dynamics') && strcmp(opts.dynamics,'lunarStatic');
assert(~static || nargin >= 9,'LuGRE:BatchLunar','lunarStatic 需要月固系输入 lunar。');
nState = 6 - 3*static;
epochs = epochs(:);
[used, signals, reference, signal] = selectSignals(obs, opts);
epochIndex = obs.epochIndex(used);
satPosition = obs.satPosition(used,:);
range = obs.range(used);
weight = obs.weight(used);
tau = epochs(epochIndex) - tRef;
halfSpan = max(max(abs(epochs(unique(epochIndex))-tRef)),1);
clockSpan = [min(tau), max(tau)];
clockBasis = clockDesign(tau, clockSpan, halfSpan, opts, static);
clockKeep = true(1,size(clockBasis,2));
if static
    clockKeep = full(any(clockBasis,1)); % 长窗口的观测空档中可能没有任何观测
    clockBasis = clockBasis(:,clockKeep);
end
signal = signal(used);
biasSignals = signals(signals ~= reference);
bias = zeros(numel(range),numel(biasSignals));
for k = 1:numel(biasSignals), bias(:,k) = signal == biasSignals(k); end
nParam = nState + size(clockBasis,2) + size(bias,2);
assert(numel(range) > nParam,'LuGRE:BatchObservations','观测数 %d 不足以估计 %d 个参数。', ...
    numel(range),nParam);

state = state0(:);
solution = struct('method',method,'converged',false,'iterations',0,'weakSearch',[],'onBoundary',false);
for outer = 1:opts.maxOuterIterations
    if static
        [A, y] = linearizeStatic(state, lunar, epochIndex, satPosition, range, clockBasis, bias);
    else
        [A, y] = linearize(state, moon0, tRef, epochs, rotations, epochIndex, satPosition, ...
            range, tau, clockBasis, bias, cfg);
    end
    scaled = sqrt(weight).*(A./columnNorm(A));
    assert(fullColumnRank(scaled, nParam),'LuGRE:BatchRank','批处理法方程秩亏，几何或参数设置不可观测。');
    if strcmp(method,'ls')
        p = solveWeighted(A, y, weight);
    else
        p = irlsL1(A, y, weight, [], opts);
        if outer <= opts.weakSearch.outerIterations
            [p, search] = weakDirectionSearch(A, y, weight, p, opts);
            solution.weakSearch = search;
            solution.onBoundary = solution.onBoundary || search.onBoundary;
        end
    end
    state = state + p(1:nState);
    solution.iterations = outer;
    if norm(p(1:3)) < opts.positionToleranceM && (static || norm(p(4:6))*halfSpan < opts.positionToleranceM)
        solution.converged = true;
        break;
    end
end
if ~solution.converged
    warning('LuGRE:BatchNotConverged','%s 批处理 %d 次外迭代未收敛。',upper(method),outer);
end
% 在最终状态重新线性化，线性参数（钟差、偏差）与最终轨道一致。
if static
    [A, y] = linearizeStatic(state, lunar, epochIndex, satPosition, range, clockBasis, bias);
else
    [A, y] = linearize(state, moon0, tRef, epochs, rotations, epochIndex, satPosition, ...
        range, tau, clockBasis, bias, cfg);
end
if strcmp(method,'ls')
    p = solveWeighted(A, y, weight);
else
    p = irlsL1(A, y, weight, p, opts);
end
state = state + p(1:nState); % 剩余轨道修正低于收敛阈值，与线性化残差一致
residual = y - A*p;
if strcmp(method,'ls')
    solution.cost = sum(weight.*residual.^2);
    scale = sqrt(solution.cost/(numel(residual)-nParam));
else
    solution.cost = sum(weight.*abs(residual));
    scale = 1.4826*median(sqrt(weight).*abs(residual)); % 稳健尺度，不受粗差支配
end
% 形式协方差（按残差尺度），列缩放后求逆以保持数值稳定。
norms = columnNorm(A);
scaled = sqrt(weight).*(A./norms);
covariance = scale^2 * (inv(full(scaled'*scaled))./(norms'*norms));

m = numel(epochs);
epochTau = epochs - tRef;
if static
    % 月固系到 ECEF 是正交变换，三维形式标准差各历元相同；速度未输出。
    positionEcef = lunar.moonEcef + permute(pagemtimes(lunar.moonFixedToEcef, state), [3 1 2]);
    velocityIcrf = nan(m,3);
    formalStd = repmat(sqrt(trace(covariance(1:3,1:3))), m, 1);
else
    [positionIcrf, velocityIcrf] = propagateArc(state, moon0, tRef, epochs, cfg);
    positionEcef = zeros(m,3); formalStd = zeros(m,1);
    for k = 1:m
        positionEcef(k,:) = (rotations(:,:,k)*positionIcrf(k,:)')';
        map = rotations(:,:,k)*[eye(3), epochTau(k)*eye(3)];
        formalStd(k) = sqrt(trace(map*covariance(1:6,1:6)*map'));
    end
end
clockColumns = nState+1:nState+size(clockBasis,2);
solution.dynamics = ternary(static,'lunarStatic','freeFlight');
solution.state = state;
solution.tRef = tRef;
solution.positionEcef = positionEcef;
solution.velocityIcrf = velocityIcrf;
epochBasis = clockDesign(epochTau, clockSpan, halfSpan, opts, static);
solution.clockM = full(epochBasis(:,clockKeep)*p(clockColumns));
solution.clockM(epochTau < clockSpan(1) | epochTau > clockSpan(2)) = NaN; % 不外推钟差
solution.clockM(full(any(epochBasis(:,~clockKeep),2))) = NaN; % 依赖未估计样条列的历元（观测空档）
solution.clockCoefficientsM = p(clockColumns);
solution.referenceSignal = reference;
solution.biasSignals = biasSignals(:);
solution.signalBiasM = p(clockColumns(end)+1:end);
solution.formalStd3dM = formalStd;
solution.robustSigmaM = scale;
solution.used = used;
solution.residualM = nan(numel(obs.range),1);
solution.residualM(used) = residual;
solution.observationCount = numel(residual);
solution.parameterCount = nParam;
end

function basis = clockDesign(tau, span, halfSpan, opts, sparseOutput)
% 钟差基函数。'polynomial'：(tau/halfSpan)^0..clockDegree；
% 'spline'：均匀节点三次 B 样条（二阶导连续），节点间距 clockKnotSpacingS。
% sparseOutput=true 时返回稀疏矩阵（样条每行只有 4 个非零元，直接按非零元构造）。
tau = tau(:);
switch opts.clockModel
    case 'polynomial'
        basis = (tau/halfSpan).^(0:opts.clockDegree);
        if sparseOutput, basis = sparse(basis); end
    case 'spline'
        spacing = opts.clockKnotSpacingS;
        % ceil+3 保证每个基函数在 [span(1), span(2)] 内都有非零支撑。
        count = ceil((span(2)-span(1))/spacing) + 3;
        if sparseOutput
            basis = sparseSpline(tau, span(1), spacing, count);
            return;
        end
        u = (tau - span(1))/spacing + 3 - (0:count-1); % 第 i 个基函数的局部坐标
        basis = zeros(size(u));
        piece = floor(u);
        a = piece==0; basis(a) = u(a).^3/6;
        b = piece==1; basis(b) = (-3*u(b).^3 + 12*u(b).^2 - 12*u(b) + 4)/6;
        c = piece==2; basis(c) = (3*u(c).^3 - 24*u(c).^2 + 60*u(c) - 44)/6;
        d = piece==3; basis(d) = (4-u(d)).^3/6;
    otherwise
        error('LuGRE:Config','clockModel 必须为 polynomial 或 spline。');
end
end

function basis = sparseSpline(tau, start, spacing, count)
% 与 clockDesign 的稠密样条逐元相同：第 i 个基函数局部坐标 u = (tau-start)/spacing + 3 - i，
% 非零的只有 floor(u) 为 0..3 的 4 个 i。
shifted = (tau - start)/spacing + 3;
column = floor(shifted) - (0:3);           % 0 起的基函数编号
u = shifted - column;
piece = floor(u);
value = zeros(size(u));
a = piece==0; value(a) = u(a).^3/6;
b = piece==1; value(b) = (-3*u(b).^3 + 12*u(b).^2 - 12*u(b) + 4)/6;
c = piece==2; value(c) = (3*u(c).^3 - 24*u(c).^2 + 60*u(c) - 44)/6;
d = piece==3; value(d) = (4-u(d)).^3/6;
row = repmat((1:numel(tau))', 1, 4);
keep = column >= 0 & column < count & piece >= 0 & piece <= 3;
basis = sparse(row(keep), column(keep)+1, value(keep), numel(tau), count);
end

function [used, signals, reference, signal] = selectSignals(obs, opts)
% 观测太少的信号无法可靠估计偏差，直接不参与，并给出数量。无身份信息记为 -1。
used = true(size(obs.range));
signal = obs.signalId(:);
signal(isnan(signal)) = -1;
signals = unique(signal);
if opts.signalBias
    counts = arrayfun(@(s)nnz(signal==s),signals);
    few = signals(counts < opts.minSignalObservations);
    if ~isempty(few)
        used = ~ismember(signal,few);
        warning('LuGRE:BatchSignal','信号 %s 观测少于 %d 条，共剔除 %d 条。', ...
            mat2str(few'),opts.minSignalObservations,nnz(~used));
    end
    signals = unique(signal(used));
    [~,most] = max(arrayfun(@(s)nnz(signal(used)==s),signals));
    reference = signals(most);
else
    reference = NaN; signals = zeros(0,1); % 不估计信号偏差
end
end

function [A, y] = linearize(state, moon0, tRef, epochs, rotations, epochIndex, satPosition, ...
    range, tau, clockBasis, bias, cfg)
% 轨道部分在当前状态处线性化：r(t) ≈ r_nom(t) + dr0 + dv0*tau。
% 对数十分钟弧段，引力梯度项引起的偏导误差约 1e-5 量级，由外迭代消除。
positionIcrf = propagateArc(state, moon0, tRef, epochs, cfg);
m = numel(epochs);
positionEcef = zeros(m,3);
for k = 1:m, positionEcef(k,:) = (rotations(:,:,k)*positionIcrf(k,:)')'; end
los = positionEcef(epochIndex,:) - satPosition;
distance = vecnorm(los,2,2);
unit = los./distance;
g = zeros(size(unit));
for k = unique(epochIndex)'
    rows = epochIndex == k;
    g(rows,:) = unit(rows,:)*rotations(:,:,k); % d(range)/d(r_ICRF)
end
A = [g, g.*tau, clockBasis, bias];
y = range - distance;
end

function [A, y] = linearizeStatic(p, lunar, epochIndex, satPosition, range, clockBasis, bias)
% 月固系静止模型：x_k = c_k + M_k*p 对 p 线性，偏导 d(range)/dp = 视线单位向量 * M_k，不积分、无 STM。
positionEcef = lunar.moonEcef + permute(pagemtimes(lunar.moonFixedToEcef, p), [3 1 2]);
los = positionEcef(epochIndex,:) - satPosition;
distance = vecnorm(los,2,2);
unit = los./distance;
g = permute(pagemtimes(permute(unit,[3 2 1]), lunar.moonFixedToEcef(:,:,epochIndex)), [3 2 1]);
A = [sparse(g), clockBasis, sparse(bias)];
y = range - distance;
end

function ok = fullColumnRank(scaled, nParam)
% 稠密：rank（SVD）；稀疏（lunarStatic 长窗口）：QR 的 R 对角元，阈值与 rank 的默认阈值同量级。
if issparse(scaled)
    R = qr(scaled, 0);
    d = abs(full(diag(R)));
    ok = numel(d) == nParam && min(d) > max(size(scaled))*eps(max(d));
else
    ok = rank(scaled) == nParam;
end
end

function p = solveWeighted(A, y, weight)
% 列缩放后用 QR（反斜杠）求解，避免形成病态正规方程。
scale = columnNorm(A);
s = sqrt(weight);
p = ((s.*(A./scale)) \ (s.*y)) ./ scale';
end

function p = irlsL1(A, y, weight, p, opts)
% 迭代重加权最小二乘求 min sum(w|y-Ap|)；epsilon 防止零残差除零。
if isempty(p), p = solveWeighted(A, y, weight); end
for iteration = 1:opts.irlsIterations
    residual = y - A*p;
    next = solveWeighted(A, y, weight./max(abs(residual),opts.irlsEpsilonM));
    change = max(abs(A*(next-p)));
    p = next;
    if change < opts.irlsToleranceM, break; end
end
end

function [p, search] = weakDirectionSearch(A, y, weight, p, opts)
% 远距离几何下径向位置与钟差几乎不可分，代价沿该方向很平且可能多峰。
% 沿法方程最弱特征方向做一维网格（每点对其余方向重新做 L1 最优），取全局最小后再局部精化。
scale = columnNorm(A);
residual = y - A*p;
omega = weight./max(abs(residual),opts.irlsEpsilonM);
scaled = A./scale;
normal = full(scaled'*(omega.*scaled));
[vectors, values] = eig((normal+normal')/2);
[~,weakest] = min(diag(values));
direction = vectors(:,weakest)./scale';
direction = direction/norm(direction(1:3)); % 位置分量为 1 m 时的参数变化
complement = null(direction');
offsets = -opts.weakSearch.halfWidthM:opts.weakSearch.stepM:opts.weakSearch.halfWidthM;
cost = nan(size(offsets));
coefficients = zeros(size(complement,2),numel(offsets));
reduced = A*complement;
searchOpts = opts; searchOpts.irlsIterations = min(opts.irlsIterations,30);
z = zeros(size(complement,2),1);
for j = 1:numel(offsets)
    shifted = y - A*(p + offsets(j)*direction);
    z = irlsL1(reduced, shifted, weight, z, searchOpts); % 相邻搜索点热启动
    cost(j) = sum(weight.*abs(shifted - reduced*z));
    coefficients(:,j) = z;
end
[~,best] = min(cost);
p = irlsL1(A, y, weight, p + offsets(best)*direction + complement*coefficients(:,best), opts);
search = struct('offsetsM',offsets,'cost',cost,'bestOffsetM',offsets(best), ...
    'direction',direction,'onBoundary',best==1 || best==numel(offsets));
end

function n = columnNorm(A)
if issparse(A)
    n = full(sqrt(sum(A.^2,1)));
else
    n = vecnorm(A,2,1);
end
n(n==0) = 1;
end

function out = ternary(condition, a, b)
if condition, out = a; else, out = b; end
end
