function positionIcrf = moonFixedToIcrf(positionMoonFixed, moonPositionIcrf, rotation)
% MOONFIXEDTOICRF 月心月固系位置转为地心 ICRF：r = r_moon + R'*p（单位 m）。
% positionMoonFixed：n x 3，或 1 x 3（各时刻同一固定点）；moonPositionIcrf：n x 3；
% rotation：icrfToMoonFixedRotation 的 3x3xn（ICRF -> 月固系）。
n = size(rotation,3);
fixed = permute(positionMoonFixed, [2 3 1]);
if size(fixed,3) == 1, fixed = repmat(fixed, 1, 1, n); end
positionIcrf = moonPositionIcrf + permute(pagemtimes(rotation, 'transpose', fixed, 'none'), [3 1 2]);
end
