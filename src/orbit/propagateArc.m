function [positionIcrf, velocityIcrf, moonState] = propagateArc(state0, moon0, tRef, times, cfg)
% PROPAGATEARC 由 tRef 的地心 ICRF 状态 [r;v] 向前、向后积分到 times（GPST s）。
% moon0：tRef 时刻地心 ICRF 月球 [r;v]（m, m/s）。times 须互不相同。
% moonState：同时积分得到的月球 [r v]（n×6）。
times = times(:);
positionIcrf = nan(numel(times),3); velocityIcrf = nan(numel(times),3);
moonState = nan(numel(times),6);
absTol = repmat([repmat(cfg.orbit.positionAbsToleranceM,3,1); ...
    repmat(cfg.orbit.velocityAbsToleranceMps,3,1)],2,1);
options = odeset('RelTol',cfg.orbit.relativeTolerance,'AbsTol',absTol);
x0 = [state0(:); moon0(:)];
for direction = [-1 1]
    index = find(sign(times-tRef) == direction);
    if isempty(index), continue; end
    [~,order] = sort(abs(times(index)-tRef));
    index = index(order);
    span = [0; times(index)-tRef];
    [~,x] = ode45(@(t,x)orbitDynamics(t,x,tRef,cfg), span, x0, options);
    if numel(span) == 2, x = x([1 end],:); end % 两点时 ode45 返回全部积分步
    positionIcrf(index,:) = x(2:end,1:3);
    velocityIcrf(index,:) = x(2:end,4:6);
    moonState(index,:) = x(2:end,7:12);
end
atRef = times == tRef;
moonState(atRef,:) = repmat(moon0(:)',nnz(atRef),1);
positionIcrf(atRef,:) = repmat(state0(1:3)',nnz(atRef),1);
velocityIcrf(atRef,:) = repmat(state0(4:6)',nnz(atRef),1);
end
