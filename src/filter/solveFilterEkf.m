function solution = solveFilterEkf(meas, epochs, rotations, rates, moon, method, cfg)
% SOLVEFILTEREKF 伪距 + Doppler 扩展卡尔曼滤波，移植自 EKF-LuGRE（main5.m 与 main5_cv_orbit_weighted_jpl.m）。
% 状态 [r; v; 钟差 c·δt (m); 钟漂 (m/s)]，8 维。
%   ekfOrbit  ：地心 ICRF 中积分地球 + 月球（可选 J2）动力学及状态转移矩阵做预测，再做伪距/Doppler 量测更新。
%   ekfCvBlend：ECEF 匀速模型 EKF 更新后，与“上一历元加权解经轨道动力学外推”的预测按 kalmanWeight 直接加权；
%               钟差状态取 EKF；协方差 = EKF 协方差 + w(1-w)·dd'（两者不一致时不至过度自信）。
% meas：每历元 buildFilterMeasurements 输出（cell）；rotations/rates：ICRF->ITRF 旋转及其导数（3x3xm）；
% moon：struct gpsSeconds（列）与 state（n×6，地心 ICRF 月球 [r v]，m、m/s）。
% 初值：首个卫星数 >= minInitSatellites 的历元做伪距 LS（自地心起算）+ Doppler LS，不用 NAV；
% 初始协方差 initCovariance='fixed' 为原项目的固定对角阵，'ls' 取 LS 形式协方差。
p = cfg.filter.(method);
m = numel(epochs);
solution = struct('method',method,'positionEcef',nan(m,3),'velocityEcef',nan(m,3), ...
    'clockM',nan(m,1),'clockDriftMps',nan(m,1),'formalStd3dM',nan(m,1), ...
    'satCount',cellfun(@(o)numel(o.range),meas),'acceptedCount',zeros(m,1), ...
    'rangeInnovationRmsM',nan(m,1),'rangeRateInnovationRmsMps',nan(m,1),'initialEpoch',NaN);
if strcmp(method,'ekfCvBlend')
    solution.kalmanPositionEcef = nan(m,3);
    solution.orbitPredictionEcef = nan(m,3);
end
sigmaOf = @(o) addSigma(o, p, cfg.filter);
k0 = find(solution.satCount >= cfg.filter.minInitSatellites, 1);
assert(~isempty(k0),'LuGRE:FilterInit','没有卫星数 >= %d 的历元，无法初始化滤波。',cfg.filter.minInitSatellites);
[x, P] = initializeState(sigmaOf(meas{k0}), cfg.filter);
solution.initialEpoch = k0;
solution.acceptedCount(k0) = 2*solution.satCount(k0);
inertial = strcmp(method,'ekfOrbit');
if inertial
    [x, P] = ecefToInertialFilter(x, P, rotations(:,:,k0), rates(:,:,k0));
end
solution = store(solution, k0, x, P, inertial, rotations(:,:,k0), rates(:,:,k0));
for k = k0+1:m
    dt = epochs(k) - epochs(k-1);
    o = sigmaOf(meas{k});
    R = rotations(:,:,k); Rd = rates(:,:,k);
    if inertial
        [x(1:6), Phi] = propagateOrbit(x(1:6), moonAt(moon, epochs(k-1)), epochs(k-1), dt, cfg, true);
        F = eye(8); F(1:6,1:6) = Phi; F(7,8) = dt;
        x(7) = x(7) + dt*x(8);
        P = symmetric(F*P*F' + processNoise(dt, p));
        if ~isempty(o.range)
            % 卫星状态转到 ICRF：r_i = R' r_e，v_i = R'(v_e - Ṙ r_i)
            o.satPosition = o.satPosition*R;
            o.satVelocity = (o.satVelocity - o.satPosition*Rd')*R;
        end
        [x, P, info] = measurementUpdate(x, P, o, cfg.filter.gateSigma);
    else
        previous = x;
        F = eye(8); F(1:3,4:6) = eye(3)*dt; F(7,8) = dt;
        [kalman, P, info] = measurementUpdate(F*x, symmetric(F*P*F' + processNoise(dt, p)), o, cfg.filter.gateSigma);
        % 轨道分支：上一历元加权解 -> ICRF -> 动力学外推 dt -> ECEF
        R0 = rotations(:,:,k-1);
        r0 = R0'*previous(1:3);
        v0 = R0'*(previous(4:6) - rates(:,:,k-1)*r0);
        propagated = propagateOrbit([r0; v0], moonAt(moon, epochs(k-1)), epochs(k-1), dt, cfg, false);
        orbit = [R*propagated(1:3); R*propagated(4:6) + Rd*propagated(1:3)];
        w = p.kalmanWeight;
        x = kalman;
        x(1:6) = w*kalman(1:6) + (1-w)*orbit;
        d = kalman(1:6) - orbit;
        P(1:6,1:6) = P(1:6,1:6) + w*(1-w)*(d*d');
        P = symmetric(P);
        solution.kalmanPositionEcef(k,:) = kalman(1:3)';
        solution.orbitPredictionEcef(k,:) = orbit(1:3)';
    end
    solution.acceptedCount(k) = info.accepted;
    solution.rangeInnovationRmsM(k) = info.rangeRms;
    solution.rangeRateInnovationRmsMps(k) = info.rangeRateRms;
    solution = store(solution, k, x, P, inertial, R, Rd);
end
end

function o = addSigma(o, p, f)
% C/N0 越低标准差越大；C/N0 缺失取上限。
scale = 10.^((f.cn0ReferenceDbHz - o.cn0)/20);
o.rangeSigma = min(max(p.rangeSigmaM*scale, f.rangeSigmaLimitsM(1)), f.rangeSigmaLimitsM(2));
o.rangeRateSigma = min(max(p.rangeRateSigmaMps*scale, f.rangeRateSigmaLimitsMps(1)), f.rangeRateSigmaLimitsMps(2));
o.rangeSigma(~isfinite(o.rangeSigma)) = f.rangeSigmaLimitsM(2);
o.rangeRateSigma(~isfinite(o.rangeRateSigma)) = f.rangeRateSigmaLimitsMps(2);
end

function [x, P] = initializeState(o, f)
% 加权伪距 LS（自地心、中位数钟差起算）求位置与钟差，再用 Doppler LS 求速度与钟漂。
sat = o.satPosition;
pc = [0; 0; 0; median(o.range - vecnorm(sat,2,2))];
w = 1./o.rangeSigma.^2;
for iteration = 1:f.initIterations
    los = pc(1:3)' - sat;
    rho = vecnorm(los,2,2);
    H = [los./rho, ones(numel(rho),1)];
    step = (H'*(w.*H)) \ (H'*(w.*(o.range - rho - pc(4))));
    pc = pc + step;
    if norm(step) < f.initStepToleranceM, break; end
end
u = (sat - pc(1:3)')./vecnorm(sat - pc(1:3)',2,2);
good = isfinite(o.rangeRate);
assert(nnz(good) >= 4,'LuGRE:FilterInit','初始化历元有效 Doppler 不足 4 个。');
H = [-u(good,:), ones(nnz(good),1)];
w = 1./o.rangeRateSigma(good).^2;
vd = (H'*(w.*H)) \ (H'*(w.*(o.rangeRate(good) - sum(u(good,:).*o.satVelocity(good,:),2))));
x = [pc(1:3); vd(1:3); pc(4); vd(4)];
if strcmp(f.initCovariance,'ls')
    % LS 形式协方差：保留径向-钟差（及径向速度-钟漂）的强相关，远距离时不会过度自信。
    Hp = [-(sat - pc(1:3)')./vecnorm(sat - pc(1:3)',2,2), ones(size(sat,1),1)];
    Cp = inv(Hp'*(Hp./o.rangeSigma.^2));
    Cv = inv(H'*(w.*H));
    P = zeros(8);
    P([1:3 7],[1:3 7]) = Cp;
    P([4:6 8],[4:6 8]) = Cv;
    P = symmetric(P);
else
    P = diag(f.initStd(:).^2);
end
end

function [x, P] = ecefToInertialFilter(x, P, R, Rd)
% r_i = R' r_e；v_i = R' v_e - R' Ṙ R' r_e
T = eye(8);
T(1:3,1:3) = R'; T(4:6,4:6) = R'; T(4:6,1:3) = -R'*Rd*R';
x = T*x; P = symmetric(T*P*T');
end

function solution = store(solution, k, x, P, inertial, R, Rd)
if inertial
    solution.positionEcef(k,:) = (R*x(1:3))';
    solution.velocityEcef(k,:) = (R*x(4:6) + Rd*x(1:3))';
else
    solution.positionEcef(k,:) = x(1:3)';
    solution.velocityEcef(k,:) = x(4:6)';
end
solution.clockM(k) = x(7);
solution.clockDriftMps(k) = x(8);
solution.formalStd3dM(k) = sqrt(trace(P(1:3,1:3))); % 迹不随旋转变化
end

function Q = processNoise(dt, p)
% 各轴白噪声加速度（PSD p.accelerationPsd）与钟漂随机游走（PSD p.clockDriftPsd）。
block = [dt^3/3, dt^2/2; dt^2/2, dt];
Q = zeros(8);
for axis = 1:3
    Q([axis axis+3],[axis axis+3]) = p.accelerationPsd*block;
end
Q([7 8],[7 8]) = p.clockDriftPsd*block;
end

function [x, P, info] = measurementUpdate(x, P, o, gateSigma)
% 伪距 + 变化率联合更新；归一化新息超过 gateSigma 的观测剔除；Joseph 形式更新协方差。
info = struct('accepted',0,'rangeRms',NaN,'rangeRateRms',NaN);
n = numel(o.range);
if n == 0, return; end
toSat = o.satPosition - x(1:3)';
rho = vecnorm(toSat,2,2);
u = toSat./rho;
relative = o.satVelocity - x(4:6)';
Hr = zeros(n,8); Hr(:,1:3) = -u; Hr(:,7) = 1;
Hd = zeros(n,8); Hd(:,4:6) = -u; Hd(:,8) = 1;
% ∂(u·v_rel)/∂r = -v_rel'(I - u u')/ρ
Hd(:,1:3) = -(relative - sum(relative.*u,2).*u)./rho;
H = [Hr; Hd];
innovation = [o.range - rho - x(7); o.rangeRate - sum(u.*relative,2) - x(8)];
variance = [o.rangeSigma.^2; o.rangeRateSigma.^2];
normalized = abs(innovation)./sqrt(max(sum((H*P).*H,2) + variance, eps));
accepted = isfinite(normalized) & normalized <= gateSigma;
info.accepted = nnz(accepted);
info.rangeRms = rmsFinite(innovation(accepted(1:n)));
info.rangeRateRms = rmsFinite(innovation([false(n,1); accepted(n+1:end)]));
if ~any(accepted), return; end
H = H(accepted,:); Rm = diag(variance(accepted));
K = (P*H') / (H*P*H' + Rm);
x = x + K*innovation(accepted);
A = eye(8) - K*H;
P = symmetric(A*P*A' + K*Rm*K');
end

function value = rmsFinite(values)
if isempty(values), value = NaN; else, value = sqrt(mean(values.^2)); end
end

function [state, Phi] = propagateOrbit(state, moonState, t0, dt, cfg, withStm)
% 定步长 RK4（子步不超过 cfg.filter.maxStepS）积分航天器与月球 [r;v]，动力学同 orbitDynamics；
% withStm 时同时积分 6x6 状态转移矩阵（梯度只含地球、月球点质量项）。
steps = max(1, ceil(abs(dt)/cfg.filter.maxStepS));
h = dt/steps;
y = [state(:); moonState(:)];
if withStm, y = [y; reshape(eye(6),36,1)]; end
f = @(t,y) derivative(t, y, t0, cfg, withStm);
t = 0;
for step = 1:steps
    k1 = f(t, y); k2 = f(t+h/2, y+h/2*k1); k3 = f(t+h/2, y+h/2*k2); k4 = f(t+h, y+h*k3);
    y = y + h*(k1 + 2*k2 + 2*k3 + k4)/6;
    t = t + h;
end
state = y(1:6);
Phi = [];
if withStm, Phi = reshape(y(13:end),6,6); end
end

function dy = derivative(t, y, t0, cfg, withStm)
dy = orbitDynamics(t, y(1:12), t0, cfg);
if ~withStm, return; end
r = y(1:3); relative = r - y(7:9);
muE = cfg.constants.muEarth; muM = cfg.constants.muMoon;
G = muE*(3*(r*r')/norm(r)^5 - eye(3)/norm(r)^3) ...
    + muM*(3*(relative*relative')/norm(relative)^5 - eye(3)/norm(relative)^3);
A = [zeros(3), eye(3); G, zeros(3)];
dy = [dy; reshape(A*reshape(y(13:end),6,6),36,1)];
end

function state = moonAt(moon, t)
state = interp1(moon.gpsSeconds, moon.state, t, 'pchip', NaN)';
assert(all(isfinite(state)),'LuGRE:MoonCoverage','月球星历未覆盖 GPS 时间 %.3f。',t);
end

function P = symmetric(P)
P = (P + P')/2;
end
