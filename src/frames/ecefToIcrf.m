function [position, velocity] = ecefToIcrf(positionEcef, velocityEcef, gpsSeconds, frame)
% ECEFTOICRF ITRF 位置(m)、速度(m/s) 转换为地心 ICRF/GCRF。
[rotation, rate] = icrfToItrfRotation(gpsSeconds, frame);
position = rotation' * positionEcef(:);
velocity = rotation' * (velocityEcef(:) - rate * position);
end
