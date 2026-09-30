function [position, velocity, clockM] = computeSatellitePvt(eph, transmitGps, mu, omegaEarth, c)
% COMPUTESATELLITEPVT 由广播星历计算卫星 ECEF 位置 (m)、速度 (m/s) 与钟差 c*dts (m)。
% eph：readRinexNav 的一条记录（Toc/Toe 用 GPST 总秒数）；transmitGps：信号发射时刻 GPST (s)。
% mu、omegaEarth、c：该星座的引力常数、地球自转角速度与光速。
sqrtA = eph.SqrtA;
e = eph.e;
tk = transmitGps - eph.ToeTotalSeconds;

% 平近点角与开普勒方程
A = sqrtA^2;
n = sqrt(mu / A^3) + eph.DeltaN;
Mk = eph.M0 + n * tk;
Ek = Mk;
for iteration = 1:20
    previous = Ek;
    Ek = Mk + e * sin(Ek);
    if abs(Ek - previous) < 1e-14, break; end
end

% 卫星钟差（含相对论项），以 Toc 为参考
F = -2 * sqrt(mu) / c^2;
relativistic = F * e * sqrtA * sin(Ek);
tc = transmitGps - eph.TocTotalSeconds;
clockM = (eph.ClockBias + eph.ClockDrift * tc + eph.ClockDriftRate * tc^2 + relativistic)*c;

% 位置
sinV = (sqrt(1 - e^2) * sin(Ek)) / (1 - e * cos(Ek));
cosV = (cos(Ek) - e) / (1 - e * cos(Ek));
vk = atan2(sinV, cosV);
phi = vk + eph.omega;
sin2Phi = sin(2 * phi);
cos2Phi = cos(2 * phi);
du = eph.Cus * sin2Phi + eph.Cuc * cos2Phi;
dr = eph.Crs * sin2Phi + eph.Crc * cos2Phi;
di = eph.Cis * sin2Phi + eph.Cic * cos2Phi;
uk = phi + du;
rk = A * (1 - e * cos(Ek)) + dr;
ik = eph.i0 + eph.IDOT * tk + di;
xPlane = rk * cos(uk);
yPlane = rk * sin(uk);
OmegaK = eph.Omega0 + (eph.OmegaDot - omegaEarth) * tk - omegaEarth * eph.Toe;
position = [xPlane * cos(OmegaK) - yPlane * cos(ik) * sin(OmegaK);
    xPlane * sin(OmegaK) + yPlane * cos(ik) * cos(OmegaK);
    yPlane * sin(ik)];

% 速度
EkDot = n / (1 - e*cos(Ek));
vkDot = sqrt(1-e^2)*EkDot/(1-e*cos(Ek));
duDot = 2*(eph.Cus*cos2Phi - eph.Cuc*sin2Phi)*vkDot;
drDot = 2*(eph.Crs*cos2Phi - eph.Crc*sin2Phi)*vkDot;
diDot = 2*(eph.Cis*cos2Phi - eph.Cic*sin2Phi)*vkDot;
uDot = vkDot + duDot;
rDot = A*e*sin(Ek)*EkDot + drDot;
iDot = eph.IDOT + diDot;
xPlaneDot = rDot*cos(uk) - rk*sin(uk)*uDot;
yPlaneDot = rDot*sin(uk) + rk*cos(uk)*uDot;
OmegaDotK = eph.OmegaDot - omegaEarth;
velocity = [xPlaneDot*cos(OmegaK) ...
    - yPlaneDot*cos(ik)*sin(OmegaK) ...
    + yPlane*sin(ik)*sin(OmegaK)*iDot ...
    - (xPlane*sin(OmegaK) + yPlane*cos(ik)*cos(OmegaK))*OmegaDotK;
    xPlaneDot*sin(OmegaK) ...
    + yPlaneDot*cos(ik)*cos(OmegaK) ...
    - yPlane*sin(ik)*cos(OmegaK)*iDot ...
    + (xPlane*cos(OmegaK) - yPlane*cos(ik)*sin(OmegaK))*OmegaDotK;
    yPlaneDot*sin(ik) + yPlane*cos(ik)*iDot];
end
