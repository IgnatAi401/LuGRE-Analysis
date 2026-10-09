function meas = buildFilterMeasurements(data, cfg)
% BUILDFILTERMEASUREMENTS 卡尔曼滤波一个历元的伪距 + Doppler 观测（ECEF，接收时刻地固系）。
% data：selectEpochMeasurements 的输出（1:3 卫星位置，4:6 卫星速度，7 修正伪距，8 C/N0，9 Doppler，14 signalId）。
% 变化率 = -λ·Doppler；Sagnac 开启时卫星位置与速度一起旋转到接收时刻地固系。
% 标准差按 C/N0 缩放：σ = σ_ref·10^((ref-C/N0)/20)，再限幅；σ_ref 由各方法自己给（见 solveFilterEkf）。
meas.satPosition = zeros(0,3); meas.satVelocity = zeros(0,3);
meas.range = zeros(0,1); meas.rangeRate = zeros(0,1);
meas.cn0 = zeros(0,1); meas.signalId = zeros(0,1);
if isempty(data), return; end
valid = all(isfinite(data(:,[1:3,7])),2) & data(:,7) > 0;
data = data(valid,:);
position = data(:,1:3); velocity = data(:,4:6);
if cfg.enable.sagnac
    angle = data(:,7) / cfg.constants.c * cfg.constants.omegaEarth;
    c = cos(angle); s = sin(angle);
    position = [c.*position(:,1)+s.*position(:,2), -s.*position(:,1)+c.*position(:,2), position(:,3)];
    velocity = [c.*velocity(:,1)+s.*velocity(:,2), -s.*velocity(:,1)+c.*velocity(:,2), velocity(:,3)];
end
meas.satPosition = position;
meas.satVelocity = velocity;
meas.range = data(:,7);
meas.signalId = data(:,14);
meas.cn0 = data(:,8);
% GPS L1 / L5，Galileo E1 / E5a / E5b
frequency = [1575.42e6; 1176.45e6; 1575.42e6; 1176.45e6; 1207.14e6];
meas.rangeRate = -cfg.constants.c ./ frequency(meas.signalId+1) .* data(:,9);
meas.rangeRate(~all(isfinite(velocity),2)) = NaN; % 无卫星速度时只用伪距
end
