function lunar = computeMoonFixedFrame(gpsSeconds, rotations, cfg, rates)
% COMPUTEMOONFIXEDFRAME 各历元的月心与月固系姿态，供月面静止位置模型使用。
% 月固系点 p 的 ECEF 位置：x_k = moonEcef(k,:)' + moonFixedToEcef(:,:,k)*p，对 p 线性。
% rotations：同一时刻 ICRF -> ITRF 的 3x3xN（icrfToItrfRotation）。月心来自 JPL Horizons（联网）。
% 给出 rates（ICRF -> ITRF 旋转的时间导数，3x3xN）时另输出 ECEF 速度项（滤波 Doppler 用）：
% v_k = moonEcefVelocity(k,:)' + moonFixedToEcefRate(:,:,k)*p；月球姿态导数用 ±derivativeStepS 中心差分。
gpsSeconds = gpsSeconds(:);
moon = fetchMoonEphemeris(gpsSeconds, cfg);
lunar.moonPositionIcrf = moon(:,2:4)*1000; % JD 列不转换
lunar.icrfToMoonFixed = icrfToMoonFixedRotation(gpsSeconds);
lunar.moonEcef = permute(pagemtimes(rotations, permute(lunar.moonPositionIcrf, [2 3 1])), [3 1 2]);
lunar.moonFixedToEcef = pagemtimes(rotations, 'none', lunar.icrfToMoonFixed, 'transpose');
if nargin < 4, return; end
moonVelocityIcrf = permute(moon(:,5:7)*1000, [2 3 1]);
lunar.moonEcefVelocity = permute(pagemtimes(rotations, moonVelocityIcrf) + ...
    pagemtimes(rates, permute(lunar.moonPositionIcrf, [2 3 1])), [3 1 2]);
h = cfg.frame.derivativeStepS;
moonRate = (icrfToMoonFixedRotation(gpsSeconds + h) - icrfToMoonFixedRotation(gpsSeconds - h))/(2*h);
lunar.moonFixedToEcefRate = pagemtimes(rates, 'none', lunar.icrfToMoonFixed, 'transpose') + ...
    pagemtimes(rotations, 'none', moonRate, 'transpose');
end
