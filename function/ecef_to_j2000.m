function [position, velocity] = ecef_to_j2000(rEcef, vEcef, gpsSeconds, frame)
% ECEF_TO_J2000 历史函数名；输出地心 ICRF/GCRF 位置(m)、速度(m/s)。
if nargin < 4
    cfg = main5_config();
    frame = cfg.frame;
end
[rotation, rate] = lugreFrameRotation(gpsSeconds, frame);
position = rotation' * rEcef(:);
velocity = rotation' * (vEcef(:) - rate * position);
end
