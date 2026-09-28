function solution = lugreLs(obs, initialPosition, initialClockM, cfg)
% LUGRELS QR/SVD 后端的反斜杠求解，不形成病态的正规方程。
% 与 DPE 共用修正观测；LS 不加轨道先验，作为纯观测基线。
solution = struct('position',[NaN;NaN;NaN], 'clockM',NaN, ...
    'valid',false,'iterations',0,'residualRmsM',NaN);
if numel(obs.range) < cfg.observation.minCount, return; end
state = [initialPosition(:);initialClockM];
weights = ones(size(obs.range));
if cfg.enable.cn0Ls, weights = obs.weight; end
for iter = 1:cfg.ls.maxIterations
    delta = state(1:3)' - obs.position;
    ranges = vecnorm(delta,2,2);
    if any(ranges == 0), return; end
    design = [delta./ranges, ones(numel(ranges),1)];
    residual = obs.range - ranges - state(4);
    a = sqrt(weights).*design;
    if rank(a,cfg.ls.rankTolerance) < 4, return; end
    update = a \ (sqrt(weights).*residual);
    if any(~isfinite(update)), return; end
    state = state + update;
    solution.iterations = iter;
    if norm(update(1:3)) <= cfg.ls.positionToleranceM && abs(update(4)) <= cfg.ls.clockToleranceM
        solution.position = state(1:3);
        solution.clockM = state(4);
        solution.valid = true;
        residual = obs.range - vecnorm(state(1:3)'-obs.position,2,2) - state(4);
        solution.residualRmsM = sqrt(mean(residual.^2));
        return;
    end
end
% 未收敛时明确保留 NaN，不把最后一次迭代当作有效解。
end
