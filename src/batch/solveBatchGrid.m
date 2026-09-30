function solution = solveBatchGrid(obs, epochs, rotations, state0, moon0, tRef, cfg)
% SOLVEBATCHGRID 多历元批处理网格 DPE：不用 IRLS，直接网格搜索加权 L1 代价。
% 未知数：tRef 的地心 ICRF [r0;v0]（6 维网格）与每个历元的钟差 c*dt_k（第 7 维）。
% 给定轨道时，历元 k 的 L1 最优钟差就是该历元残差的加权中位数，等价于无限细的钟差网格，
% 因此网格只铺 6 维轨道状态。审计见 experiments/grid7dAudit.m：
%   - 单一钟差参数不可行（10 s 内钟差偏离常数 50~180 m），每历元钟差信息损失仅 1.5~3%；
%   - 网格点位置用状态转移矩阵（STM）外推，50 km/50 m/s 偏移的误差 <3 mm，无需逐点积分；
%   - 网格沿信息矩阵的白化方向布置，每轴 pointsPerAxis 点、逐级步长减半。线性化后的代价
%     （消去钟差后）是凸函数，逐级细化收敛到全局最优；外迭代重新积分并线性化。
% 只有 1 颗星的历元钟差吸收全部残差，不含轨道信息，不参与搜索。
% 输入与 solveBatchIrls 相同；输出字段与之兼容（无信号偏差），search 字段记录网格过程。
opts = cfg.batch; g = opts.grid;
validateattributes(g.pointsPerAxis,{'double'},{'scalar','integer','odd','>=',3});
epochs = epochs(:);
m = numel(epochs);
epochIndex = obs.epochIndex(:);
weight = obs.weight(:);
counts = accumarray(epochIndex, 1, [m 1]);
search = counts(epochIndex) >= 2;
multi = nnz(counts >= 2);
assert(nnz(search) > 6 + multi, 'LuGRE:BatchObservations', ...
    '多星历元观测 %d 条，不足以估计 6 个轨道参数与 %d 个历元钟差。', nnz(search), multi);
layout = epochLayout(epochIndex(search), weight(search));
offsets = gridOffsets(g.pointsPerAxis);
tau = epochs - tRef;
halfSpan = max(max(abs(tau(counts > 0))), 1);

state = state0(:);
solution = struct('method','dpeGrid','converged',false,'iterations',0,'weakSearch',[],'onBoundary',false);
history = struct('levels',{},'points',{},'boundaryMoves',{},'shiftPositionM',{},'shiftVelocityMps',{});
timer = tic;
for outer = 1:opts.maxOuterIterations
    [H, y] = linearize(state, moon0, tRef, epochs, rotations, obs, cfg);
    T = whitening(H(search,:), layout);
    if outer == 1
        % 先验盒 |dr|<=priorPositionM、|dv|<=priorVelocityMps 映射到白化坐标的各轴半宽
        Tinv = inv(T);
        half = vecnorm(Tinv(:,1:3),2,2)*g.priorPositionM + vecnorm(Tinv(:,4:6),2,2)*g.priorVelocityMps;
    else
        half = g.refineHalfWidth*ones(6,1);
    end
    [u, info] = coarseToFine(H(search,:)*T, y(search), layout, half, offsets, g);
    dx = T*u;
    state = state + dx;
    info.shiftPositionM = norm(dx(1:3)); info.shiftVelocityMps = norm(dx(4:6));
    history(end+1) = info; %#ok<AGROW>
    solution.iterations = outer;
    solution.onBoundary = solution.onBoundary || info.levels >= g.maxLevels;
    if norm(dx(1:3)) < opts.positionToleranceM && norm(dx(4:6))*halfSpan < opts.positionToleranceM
        solution.converged = true;
        break;
    end
end
if ~solution.converged
    warning('LuGRE:BatchNotConverged','DPEGRID 批处理 %d 次外迭代未收敛。',outer);
end
if solution.onBoundary
    warning('LuGRE:GridLevels','DPEGRID 网格达到 maxLevels=%d 仍未细化到 finalStep，可能需加大先验盒或层数。',g.maxLevels);
end

% 最终状态处：每历元钟差 = 加权中位数，残差与形式协方差（每历元钟差模型）。
[H, y, Phi, positionIcrf, velocityIcrf] = linearize(state, moon0, tRef, epochs, rotations, obs, cfg);
clock = nan(m,1);
for k = find(counts > 0)'
    rows = epochIndex == k;
    clock(k) = weightedMedian(y(rows), weight(rows));
end
residual = y - clock(epochIndex);
solution.cost = sum(weight(search).*abs(residual(search)));
scale = 1.4826*median(sqrt(weight(search)).*abs(residual(search)));
T = whitening(H(search,:), layout);
covariance = scale^2 * (T*T'); % T'*(H'WH，历元内去均值)*T = I

positionEcef = zeros(m,3); formalStd = zeros(m,1);
for k = 1:m
    positionEcef(k,:) = (rotations(:,:,k)*positionIcrf(k,:)')';
    map = rotations(:,:,k)*Phi(:,:,k);
    formalStd(k) = sqrt(trace(map*covariance*map'));
end
solution.state = state;
solution.tRef = tRef;
solution.positionEcef = positionEcef;
solution.velocityIcrf = velocityIcrf;
solution.clockM = clock;
solution.clockCoefficientsM = clock;
solution.referenceSignal = NaN;
solution.biasSignals = zeros(0,1);
solution.signalBiasM = zeros(0,1);
solution.formalStd3dM = formalStd;
solution.covariance = covariance;
solution.robustSigmaM = scale;
solution.used = true(size(obs.range));
solution.residualM = residual;
solution.observationCount = nnz(search);
solution.parameterCount = 6 + multi;
solution.search = struct('history',history,'seconds',toc(timer), ...
    'pointsPerAxis',g.pointsPerAxis,'totalPoints',sum([history.points]));
end

function [H, y, Phi, positionIcrf, velocityIcrf] = linearize(state, moon0, tRef, epochs, rotations, obs, cfg)
% 伪距对 [r0;v0] 的偏导，轨道偏导用有限差分 STM（动力学在此尺度上几乎线性）。
steps = [1e3*ones(3,1); ones(3,1)];
[positionIcrf, velocityIcrf] = propagateArc(state, moon0, tRef, epochs, cfg);
m = numel(epochs);
Phi = zeros(3,6,m);
for j = 1:6
    dx = zeros(6,1); dx(j) = steps(j);
    p = propagateArc(state + dx, moon0, tRef, epochs, cfg);
    Phi(:,j,:) = permute((p - positionIcrf)/steps(j), [2 3 1]);
end
positionEcef = zeros(m,3);
for k = 1:m, positionEcef(k,:) = (rotations(:,:,k)*positionIcrf(k,:)')'; end
idx = obs.epochIndex(:);
los = positionEcef(idx,:) - obs.satPosition;
distance = vecnorm(los,2,2);
unit = los./distance;
H = zeros(numel(idx),6);
for k = unique(idx)'
    rows = idx == k;
    H(rows,:) = unit(rows,:)*rotations(:,:,k)*Phi(:,:,k);
end
y = obs.range(:) - distance;
end

function T = whitening(H, layout)
% 每历元钟差自由：H 在历元内减去加权均值后，信息矩阵 M = Hc'*W*Hc。
% 返回 T 使 T'*M*T = I；列缩放后做特征分解以避免 m 与 m/s 混合单位的病态。
Hc = H - groupWeightedMean(H, layout);
n = vecnorm(Hc,2,1); n(n==0) = 1;
Hs = Hc./n;
Ms = Hs'*(layout.weight.*Hs);
[V, D] = eig((Ms+Ms')/2);
d = diag(D);
assert(min(d) > max(d)*1e-14, 'LuGRE:BatchRank','DPEGRID 信息矩阵秩亏，轨道状态不可观测。');
T = (V./n')./sqrt(d');
end

function g = groupWeightedMean(X, layout)
wsum = accumarray(layout.group, layout.weight);
g = zeros(size(X));
for j = 1:size(X,2)
    mu = accumarray(layout.group, layout.weight.*X(:,j))./wsum;
    g(:,j) = mu(layout.group);
end
end

function [u, info] = coarseToFine(HT, y, layout, half, offsets, g)
% 以 0 为中心，每轴 P 个点；最优点不在边界的轴步长减半，在边界的轴保持步长并移动中心。
% 所有轴步长到 finalStep 且最优点在内部时停止。
reach = (g.pointsPerAxis - 1)/2;
step = max(half/reach, g.finalStep);
u = zeros(6,1);
info = struct('levels',0,'points',0,'boundaryMoves',0,'shiftPositionM',NaN,'shiftVelocityMps',NaN);
for level = 1:g.maxLevels
    U = u + step.*offsets;
    cost = gridCost(HT, y, layout, U, g.chunkSize);
    [~, best] = min(cost);
    u = U(:,best);
    edge = abs(offsets(:,best)) == reach;
    info.levels = level; info.points = info.points + size(U,2);
    info.boundaryMoves = info.boundaryMoves + any(edge);
    if all(step <= g.finalStep) && ~any(edge), break; end
    step(~edge) = max(step(~edge)/2, g.finalStep);
end
end

function cost = gridCost(HT, y, layout, U, chunkSize)
% 每个网格点：残差 r = y - HT*u，历元内取加权中位数为钟差，代价 sum(w*|r - 钟差|)。
G = size(U,2); cost = zeros(1,G);
[mEff, K] = size(layout.order);
W = layout.paddedWeight; % mEff x K，补位为 0
for first = 1:chunkSize:G
    cols = first:min(first+chunkSize-1, G);
    R = [y - HT*U(:,cols); nan(1,numel(cols))];
    Rp = reshape(R(layout.order,:), mEff, K, numel(cols));
    if layout.equalWeights
        med = median(Rp, 2, 'omitnan');
    else
        med = paddedWeightedMedian(Rp, W);
    end
    cost(cols) = reshape(sum(W.*abs(Rp - med), [1 2], 'omitnan'), 1, []);
end
end

function med = paddedWeightedMedian(Rp, W)
% Rp: mEff x K x G（补位为 NaN，排序后在末尾）；W: mEff x K 权重（补位为 0）。
[S, order] = sort(Rp, 2);
[mEff, K, G] = size(Rp);
rows = repmat((1:mEff)', 1, K, G);
Ws = W(sub2ind([mEff K], rows, order));
Ws(isnan(S)) = 0;
cumulative = cumsum(Ws, 2);
pick = sum(cumulative < cumulative(:,end,:)/2, 2) + 1;
med = S(sub2ind(size(S), (1:mEff)'.*ones(1,1,G), pick, reshape(1:G,1,1,[]).*ones(mEff,1)));
end

function layout = epochLayout(epochIndex, weight)
% 把观测按历元排成 mEff x K 的补位矩阵（补位指向末尾的 NaN 行），供网格批量计算。
[~, ~, group] = unique(epochIndex);
counts = accumarray(group, 1);
mEff = numel(counts); K = max(counts);
[~, sorter] = sort(group);
order = zeros(mEff, K);
position = 1;
for k = 1:mEff
    order(k,1:counts(k)) = sorter(position:position+counts(k)-1)';
    position = position + counts(k);
end
pad = order == 0;
paddedWeight = zeros(mEff, K);
paddedWeight(~pad) = weight(order(~pad));
order(pad) = numel(epochIndex) + 1;
layout = struct('group',group,'weight',weight,'order',order,'paddedWeight',paddedWeight, ...
    'equalWeights',all(weight == weight(1)));
end

function offsets = gridOffsets(P)
reach = (P - 1)/2;
axes = repmat({-reach:reach}, 1, 6);
[a{1:6}] = ndgrid(axes{:});
offsets = cell2mat(cellfun(@(x)x(:)', a', 'UniformOutput', false));
end

function value = weightedMedian(x, w)
[x, order] = sort(x(:)); w = w(order);
cumulative = cumsum(w);
value = x(find(cumulative >= cumulative(end)/2, 1));
end
