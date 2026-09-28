function solution = lugreDpe(obs, center, clockM, predicted, epoch, cfg)
% LUGREDPE 分块四维网格：加权 L1 伪距残差 + 可选轨道距离先验。
solution = struct('position',[NaN;NaN;NaN], 'clockM',NaN, ...
    'valid',false,'cost',NaN,'onBoundary',false);
if numel(obs.range) < cfg.observation.minCount, return; end
if epoch < cfg.dpe.switchEpoch
    offsets = cfg.dpe.coarseOffsetsM;
else
    offsets = cfg.dpe.fineOffsetsM;
end
clockOffsets = cfg.dpe.clockOffsetsM;
side = numel(offsets);
shape = [side,side,side,numel(clockOffsets)];
count = prod(shape);
weights = ones(size(obs.range));
if cfg.enable.cn0Dpe, weights = obs.weight; end
best = Inf; bestIndex = 0;
for first = 1:cfg.dpe.chunkSize:count
    index = first:min(first+cfg.dpe.chunkSize-1,count);
    [ix,iy,iz,it] = ind2sub(shape,index);
    positions = center(:) + [offsets(ix); offsets(iy); offsets(iz)];
    clocks = clockM + clockOffsets(it);
    cost = zeros(1,numel(index));
    for j = 1:numel(obs.range)
        ranges = vecnorm(positions-obs.position(j,:)',2,1);
        cost = cost + weights(j)*abs(obs.range(j)-ranges-clocks);
    end
    if cfg.enable.orbit && all(isfinite(predicted))
        cost = cost + cfg.dpe.orbitWeight*vecnorm(positions-predicted(:),2,1);
    end
    [value,local] = min(cost);
    if value < best
        best = value; bestIndex = index(local);
        solution.position = positions(:,local);
        solution.clockM = clocks(local);
    end
end
solution.valid = isfinite(best);
solution.cost = best;
if solution.valid
    [ix,iy,iz,it] = ind2sub(shape,bestIndex);
    subs = [ix,iy,iz,it];
    solution.onBoundary = any(subs == 1 | subs == shape);
end
end
