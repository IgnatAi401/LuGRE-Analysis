function raw = addSatelliteStates(raw, gps, galileo, maxAgeS)
% ADDSATELLITESTATES 为每条 RAW 测量匹配 Toe 最近的广播星历，计算发射时刻的卫星状态。
% 新增字段：satPosition、satVelocity（ECEF, m, m/s）与 satClockM（卫星钟差 c*dts, m）。
% raw 须含 txTime（GPST 发射时刻）、signalId（0/1 GPS，2~4 Galileo）、svId。
% maxAgeS：发射时刻与 Toe 的最大允许间隔(s)，默认 4 h；超出置 NaN，
% 防止误用其他日期的星历时静默得到 km 级错误的卫星位置。
if nargin < 4, maxAgeS = 4*3600; end
c = 299792458;
omegaEarth = 7.2921151467e-5;
muGps = 3.986005e14;       % IS-GPS-200
muGalileo = 3.986004418e14; % Galileo OS SIS ICD
missingCount = 0; staleCount = 0;
% PRN 与 Toe 预先取成数组，避免每条观测都从上万条星历结构体重新拼接。
[gpsPrn, gpsToe] = ephemerisKeys(gps);
[galileoPrn, galileoToe] = ephemerisKeys(galileo);
for k = 1:numel(raw)
    raw(k).satPosition = nan(3,1);
    raw(k).satVelocity = nan(3,1);
    raw(k).satClockM = NaN;
    switch raw(k).signalId
        case {0, 1}
            ephemerides = gps; prn = gpsPrn; toe = gpsToe; mu = muGps;
        case {2, 3, 4}
            ephemerides = galileo; prn = galileoPrn; toe = galileoToe; mu = muGalileo;
        otherwise
            continue; % 未知星座
    end
    candidates = find(prn == raw(k).svId);
    if isempty(candidates)
        missingCount = missingCount + 1;
        continue;
    end
    [age, best] = min(abs(toe(candidates) - raw(k).txTime));
    if age > maxAgeS
        staleCount = staleCount + 1;
        continue;
    end
    [raw(k).satPosition, raw(k).satVelocity, raw(k).satClockM] = ...
        computeSatellitePvt(ephemerides(candidates(best)), raw(k).txTime, mu, omegaEarth, c);
end
if missingCount > 0
    warning('LuGRE:EphemerisMissing','%d 条观测找不到对应卫星的星历，已置 NaN。',missingCount);
end
if staleCount > 0
    warning('LuGRE:EphemerisStale', ...
        '%d/%d 条观测的最近星历距发射时刻超过 %.1f h，已置 NaN；请核对 RINEX 日期。', ...
        staleCount, numel(raw), maxAgeS/3600);
end
end

function [prn, toe] = ephemerisKeys(ephemerides)
if isempty(ephemerides)
    prn = zeros(1,0); toe = zeros(1,0);
else
    prn = [ephemerides.PRN]; toe = [ephemerides.ToeTotalSeconds];
end
end
