function [position, velocity] = icrfToEcef(positionIcrf, velocityIcrf, gpsSeconds, frame)
% ICRFTOECEF 地心 ICRF/GCRF 位置(m)、速度(m/s) 转换为 ITRF。
[rotation, rate] = icrfToItrfRotation(gpsSeconds, frame);
position = rotation * positionIcrf(:);
velocity = rotation * velocityIcrf(:) + rate * positionIcrf(:);
end
