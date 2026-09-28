function [position, velocity] = j2000_to_ecef(rIcrf, vIcrf, gpsSeconds, frame)
% J2000_TO_ECEF 输入地心 ICRF/GCRF；输出 ITRF 位置(m)、速度(m/s)。
if nargin < 4
    cfg = main5_config();
    frame = cfg.frame;
end
[rotation, rate] = lugreFrameRotation(gpsSeconds, frame);
position = rotation * rIcrf(:);
velocity = rotation * vIcrf(:) + rate * rIcrf(:);
end
