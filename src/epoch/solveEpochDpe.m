function solution = solveEpochDpe(obs, center, clockM, predicted, epoch, cfg)
% SOLVEEPOCHDPE 单历元网格 DPE：分块三维位置网格 + 等价钟差最小化，
% 代价为加权 L1 伪距残差 + 可选轨道距离先验。epoch 为历元序号，决定粗/细网格。
solution = struct('position',[NaN;NaN;NaN], 'clockM',NaN, ...
    'valid',false,'cost',NaN,'onBoundary',false);
if numel(obs.range) < max(4,cfg.observation.minCount), return; end
if any(~isfinite(obs.satPosition),'all') || any(~isfinite(obs.range)), return; end
if cfg.enable.orbit && cfg.dpe.usePredictionCenter && all(isfinite(predicted))
    center = predicted(:);
end
delta = center(:)' - obs.satPosition;
ranges = vecnorm(delta,2,2);
if any(ranges==0) || rank([delta./ranges,ones(size(ranges))],cfg.dpe.rankTolerance)<4
    return;
end
if epoch < cfg.dpe.switchEpoch
    offsets = cfg.dpe.coarseOffsetsM;
else
    offsets = cfg.dpe.fineOffsetsM;
end
clockOffsets = cfg.dpe.clockOffsetsM;
solution = searchGrid(obs,center,clockM,predicted,offsets,clockOffsets,cfg);
if ~solution.valid, return; end
% 每级覆盖上一级的一个步长，保留原候选，因此代价不会增大。
positionRadius = max([diff(offsets),0]);
clockRadius = max([diff(clockOffsets),0]);
for step = cfg.dpe.refineStepsM
    offsets = (-ceil(positionRadius/step):ceil(positionRadius/step))*step;
    clockOffsets = (-ceil(clockRadius/step):ceil(clockRadius/step))*step;
    refined = searchGrid(obs,solution.position,solution.clockM,predicted,offsets,clockOffsets,cfg);
    boundary = solution.onBoundary || refined.onBoundary;
    if refined.valid && refined.cost <= solution.cost, solution = refined; end
    solution.onBoundary = boundary;
    positionRadius = step; clockRadius = step;
end
end

function solution = searchGrid(obs,center,clockM,predicted,offsets,clockOffsets,cfg)
solution = struct('position',nan(3,1),'clockM',NaN,'valid',false,'cost',NaN,'onBoundary',false);
side = numel(offsets);
spatialCount = side^3;
clocks = clockM + clockOffsets;
weights = ones(size(obs.range));
if cfg.enable.cn0Dpe, weights = obs.weight; end
if any(~isfinite(weights) | weights<=0), return; end
best = Inf; bestIndex = Inf;
for first = 1:cfg.dpe.chunkSize:spatialCount
    index = first:min(first+cfg.dpe.chunkSize-1,spatialCount);
    [ix,iy,iz] = ind2sub([side,side,side],index);
    positions = center(:) + [offsets(ix);offsets(iy);offsets(iz)];
    residual = zeros(numel(obs.range),numel(index));
    for j = 1:numel(obs.range)
        residual(j,:) = obs.range(j)-vecnorm(positions-obs.satPosition(j,:)',2,1);
    end
    % 固定位置时 L1 目标的极小值位于加权中位数；离散网格只需比较两邻点。
    [sorted,order] = sort(residual,1);
    cumulative = cumsum(reshape(weights(order),size(order)),1);
    medianRow = sum(cumulative < sum(weights)/2,1)+1;
    medianClock = sorted(sub2ind(size(sorted),medianRow,1:numel(index)));
    lower = ones(size(index));
    for k = 2:numel(clocks)
        lower(medianClock>=clocks(k)) = k;
    end
    upper = min(lower+1,numel(clocks));
    lowerCost = sum(weights.*abs(residual-clocks(lower)),1);
    upperCost = sum(weights.*abs(residual-clocks(upper)),1);
    useUpper = upperCost < lowerCost;
    it = lower; it(useUpper) = upper(useUpper);
    cost = min(lowerCost,upperCost);
    if cfg.enable.orbit && all(isfinite(predicted))
        cost = cost + cfg.dpe.orbitWeight*vecnorm(positions-predicted(:),2,1);
    end
    % 保持原四维枚举的平局顺序：钟差维最后、X 维最先。
    originalIndex = index + (it-1)*spatialCount;
    value = min(cost);
    tied = find(cost==value);
    [candidateIndex,tie] = min(originalIndex(tied));
    local = tied(tie);
    if value < best || (value==best && candidateIndex<bestIndex)
        best = value; bestIndex = candidateIndex;
        solution.position = positions(:,local);
        solution.clockM = clocks(it(local));
    end
end
solution.valid = isfinite(best);
solution.cost = best;
if solution.valid
    shape = [side,side,side,numel(clocks)];
    [ix,iy,iz,it] = ind2sub(shape,bestIndex);
    subs = [ix,iy,iz,it];
    solution.onBoundary = any(subs==1 | subs==shape);
end
end
