function obs = buildObservations(data, origin, cfg)
% BUILDOBSERVATIONS 所有方法共用的观测修正（Sagnac、可选旧电离层项）与 C/N0 权重。
% data：selectEpochMeasurements 的输出（至少 8 列：1:3 卫星位置，7 修正伪距，8 C/N0）。
obs.satPosition = zeros(0,3); obs.range = zeros(0,1); obs.weight = zeros(0,1);
obs.cn0 = zeros(0,1); obs.svId = zeros(0,1); obs.signalId = zeros(0,1);
if isempty(data), return; end
valid = all(isfinite(data(:,[1:3,7])),2) & data(:,7) > 0;
if cfg.enable.cn0Dpe || cfg.enable.cn0Ls
    valid = valid & isfinite(data(:,8));
end
data = data(valid,:);
obs.satPosition = data(:,1:3);
obs.range = data(:,7);
if cfg.enable.sagnac
    angle = obs.range / cfg.constants.c * cfg.constants.omegaEarth;
    x = obs.satPosition(:,1); y = obs.satPosition(:,2);
    obs.satPosition(:,1) = cos(angle).*x + sin(angle).*y;
    obs.satPosition(:,2) = -sin(angle).*x + cos(angle).*y;
end
if cfg.enable.legacyIonosphere
    los = obs.satPosition - origin(:)';
    elevation = asind(los(:,3)./vecnorm(los,2,2));
    ion = cfg.observation.legacyZenithDelayS * cfg.constants.c .* ...
        (1 + 16*(0.53-elevation/180).^3);
    obs.range = obs.range - ion;
end
% 保留身份与 C/N0，供多历元批处理估计信号间偏差和加权；8 列输入无身份信息。
obs.cn0 = data(:,8);
if size(data,2) >= 14
    obs.svId = data(:,13); obs.signalId = data(:,14);
else
    obs.svId = nan(size(obs.range)); obs.signalId = nan(size(obs.range));
end
obs.weight = ones(size(obs.range));
if cfg.enable.cn0Dpe || cfg.enable.cn0Ls
    w = 10.^((data(:,8)-cfg.observation.cn0ReferenceDbHz)/10);
    obs.weight = min(max(w,cfg.observation.weightLimits(1)),cfg.observation.weightLimits(2));
end
end
