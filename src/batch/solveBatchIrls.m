function solution = solveBatchIrls(obs, epochs, rotations, state0, moon0, tRef, method, cfg)
% SOLVEBATCHIRLS 多历元批处理定轨：整段观测联合估计一组状态。
% 未知数：tRef 的地心 ICRF [r0;v0]、钟差 c*dt(t)（整段多项式或三次 B 样条）、
% 各信号相对参考信号的常值偏差。
% method = 'dpeIrls'：加权 L1 伪距残差代价（IRLS 求解），可选沿法方程最弱方向一维搜索（weakSearch）；
% method = 'ls' ：同一模型、同一观测的加权 L2 代价（高斯-牛顿），作为公平对照。
% obs 字段：epochIndex、satPosition(ECEF, 已含 Sagnac)、range(已加卫星钟差)、weight、signalId。
% rotations(:,:,k)：epochs(k) 的 ICRF->ECEF 矩阵。单位 m、m/s、s。
opts = cfg.batch;
assert(any(strcmp(method,{'dpeIrls','ls'})),'LuGRE:BatchMethod','method 必须为 dpeIrls 或 ls。');
epochs = epochs(:);
[used, signals, reference, signal] = selectSignals(obs, opts);
epochIndex = obs.epochIndex(used);
satPosition = obs.satPosition(used,:);
range = obs.range(used);
weight = obs.weight(used);
tau = epochs(epochIndex) - tRef;
halfSpan = max(max(abs(epochs(unique(epochIndex))-tRef)),1);
clockSpan = [min(tau), max(tau)];
clockBasis = clockDesign(tau, clockSpan, halfSpan, opts);
signal = signal(used);
biasSignals = signals(signals ~= reference);
bias = zeros(numel(range),numel(biasSignals));
for k = 1:numel(biasSignals), bias(:,k) = signal == biasSignals(k); end
nParam = 6 + size(clockBasis,2) + size(bias,2);
assert(numel(range) > nParam,'LuGRE:BatchObservations','观测数 %d 不足以估计 %d 个参数。', ...
    numel(range),nParam);

state = state0(:);
solution = struct('method',method,'converged',false,'iterations',0,'weakSearch',[],'onBoundary',false);
for outer = 1:opts.maxOuterIterations
    [A, y] = linearize(state, moon0, tRef, epochs, rotations, epochIndex, satPosition, ...
        range, tau, clockBasis, bias, cfg);
    scaled = sqrt(weight).*(A./columnNorm(A));
    assert(rank(scaled) == nParam,'LuGRE:BatchRank','批处理法方程秩亏，几何或参数设置不可观测。');
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
    state = state + p(1:6);
    solution.iterations = outer;
    if norm(p(1:3)) < opts.positionToleranceM && norm(p(4:6))*halfSpan < opts.positionToleranceM
        solution.converged = true;
        break;
    end
end
if ~solution.converged
    warning('LuGRE:BatchNotConverged','%s 批处理 %d 次外迭代未收敛。',upper(method),outer);
end
% 在最终状态重新线性化，线性参数（钟差、偏差）与最终轨道一致。
[A, y] = linearize(state, moon0, tRef, epochs, rotations, epochIndex, satPosition, ...
    range, tau, clockBasis, bias, cfg);
if strcmp(method,'ls')
    p = solveWeighted(A, y, weight);
else
    p = irlsL1(A, y, weight, p, opts);
end
state = state + p(1:6); % 剩余轨道修正低于收敛阈值，与线性化残差一致
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
covariance = scale^2 * (inv(scaled'*scaled)./(norms'*norms));

[positionIcrf, velocityIcrf] = propagateArc(state, moon0, tRef, epochs, cfg);
m = numel(epochs);
positionEcef = zeros(m,3); formalStd = zeros(m,1);
epochTau = epochs - tRef;
for k = 1:m
    positionEcef(k,:) = (rotations(:,:,k)*positionIcrf(k,:)')';
    map = rotations(:,:,k)*[eye(3), epochTau(k)*eye(3)];
    formalStd(k) = sqrt(trace(map*covariance(1:6,1:6)*map'));
end
clockColumns = 7:6+size(clockBasis,2);
solution.state = state;
solution.tRef = tRef;
solution.positionEcef = positionEcef;
solution.velocityIcrf = velocityIcrf;
solution.clockM = clockDesign(epochTau, clockSpan, halfSpan, opts)*p(clockColumns);
solution.clockM(epochTau < clockSpan(1) | epochTau > clockSpan(2)) = NaN; % 不外推钟差
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

function basis = clockDesign(tau, span, halfSpan, opts)
% 钟差基函数。'polynomial'：(tau/halfSpan)^0..clockDegree；
% 'spline'：均匀节点三次 B 样条（二阶导连续），节点间距 clockKnotSpacingS。
tau = tau(:);
switch opts.clockModel
    case 'polynomial'
        basis = (tau/halfSpan).^(0:opts.clockDegree);
    case 'spline'
        spacing = opts.clockKnotSpacingS;
        % ceil+3 保证每个基函数在 [span(1), span(2)] 内都有非零支撑。
        count = ceil((span(2)-span(1))/spacing) + 3;
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
normal = scaled'*(omega.*scaled);
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
n = vecnorm(A,2,1);
n(n==0) = 1;
end
