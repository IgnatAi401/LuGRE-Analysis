function positionMoonFixed = icrfToMoonFixed(positionIcrf, moonPositionIcrf, rotation)
% ICRFTOMOONFIXED 地心 ICRF 位置 (n x 3, m) 转为月心月固系位置 (n x 3, m)。
% moonPositionIcrf：同一时刻地心 ICRF 月心位置 (n x 3, m)；rotation：icrfToMoonFixedRotation 的 3x3xn。
relative = permute(positionIcrf - moonPositionIcrf, [2 3 1]);
positionMoonFixed = permute(pagemtimes(rotation, relative), [3 1 2]);
end
