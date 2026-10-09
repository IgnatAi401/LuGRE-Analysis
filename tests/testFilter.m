function testFilter()
% TESTFILTER 扩展卡尔曼滤波 ekfOrbit / ekfCvBlend 的合成数据检查（离线，不需要 Aerospace Toolbox）。
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root,genpath(fullfile(root,'src')));
cfg = lugreConfig();
cfg.enable.j2 = false;
cfg.orbit.relativeTolerance = 1e-11;

% 合成场景：地心约 1e8 m 的自由飞行轨道，6 颗 GNSS 星每历元随机 3~6 颗可见，钟差线性漂移。
rng(11);
tRef = 1421300000;
moon0 = [3.8e8;0;0;0;1000;0];
state = [8e7;5e7;2e7;-800;1200;300];
epochs = tRef + (0:1200)';
m = numel(epochs);
[truthPos, truthVel, moonState] = propagateArc(state,moon0,tRef,epochs,cfg);
moon = struct('gpsSeconds',epochs,'state',moonState);
rotations = repmat(eye(3),1,1,m); rates = zeros(3,3,m); % ICRF 与 ECEF 取同一坐标系
radius = 2.656e7; n = 1.4585e-4;
planes = [0 55 110 165 220 275]*pi/180; phase = (0:5)*1.1;
clockM = @(t) -1e5 + 30*(t-tRef);
meas = cell(m,1);
for k = 1:m
    visible = randperm(6, 2 + randi(4));
    if k == 1, visible = 1:6; end
    sats = zeros(numel(visible),6);
    for j = 1:numel(visible)
        a = phase(visible(j)) + n*(epochs(k)-tRef);
        node = [cos(planes(visible(j))); sin(planes(visible(j))); 0];
        inPlane = cross([0;0;1],node)*cosd(55) + [0;0;1]*sind(55);
        sats(j,:) = [radius*(cos(a)*node + sin(a)*inPlane); radius*n*(-sin(a)*node + cos(a)*inPlane)]';
    end
    toSat = sats(:,1:3) - truthPos(k,:);
    u = toSat./vecnorm(toSat,2,2);
    meas{k} = struct('satPosition',sats(:,1:3),'satVelocity',sats(:,4:6), ...
        'range',vecnorm(toSat,2,2) + clockM(epochs(k)) + 5*randn(numel(visible),1), ...
        'rangeRate',sum(u.*(sats(:,4:6)-truthVel(k,:)),2) + 30 + 0.05*randn(numel(visible),1), ...
        'cn0',35*ones(numel(visible),1),'signalId',zeros(numel(visible),1));
end

for method = {'ekfOrbit','ekfCvBlend'}
    s = solveFilterEkf(meas, epochs, rotations, rates, moon, method{1}, cfg);
    err = vecnorm(s.positionEcef - truthPos,2,2);
    assert(s.initialEpoch == 1 && all(isfinite(err)), '%s 输出含 NaN', method{1});
    late = err(end-299:end);
    assert(median(late) < 0.2*err(1), '%s 末段误差 %.0f m 未明显小于初值误差 %.0f m', method{1}, median(late), err(1));
    fprintf('testFilter: %s 初值误差 %.0f m，末 300 s 误差中位数 %.0f m，形式 3D std %.0f m\n', ...
        method{1}, err(1), median(late), s.formalStd3dM(end));
    if strcmp(method{1},'ekfOrbit')
        % 动力学模型与合成真值一致时，误差应与形式精度相符。
        assert(median(late) < 3*median(s.formalStd3dM(end-299:end)), ...
            'ekfOrbit 末段误差 %.0f m 超过形式 3σ %.0f m', median(late), 3*median(s.formalStd3dM(end-299:end)));
        assert(abs(s.clockDriftMps(end) - 30) < 1, 'ekfOrbit 钟漂 %.2f m/s（真值 30）', s.clockDriftMps(end));
    end
end
fprintf('testFilter: all checks passed.\n');
end
