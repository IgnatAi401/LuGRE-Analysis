function [rotation, rate] = lugreFrameRotation(gpsSeconds, frame)
% LUGREFRAMEROTATION ICRF/GCRF -> ITRF，IAU 2000/2006；时间输入 GPST。
% GCRF 与地心 ICRF 轴对齐。使用 IERS 的 UT1、极移，不用 GMST 冒充 J2000。
if exist('dcmeci2ecef', 'file') ~= 2
    error('LuGRE:AerospaceRequired', '轨道模式需要 Aerospace Toolbox 的 dcmeci2ecef。');
end
rotation = atTime(gpsSeconds, frame);
if nargout > 1
    h = frame.derivativeStepS;
    rate = (atTime(gpsSeconds + h, frame) - atTime(gpsSeconds - h, frame)) / (2*h);
end
end

function rotation = atTime(gpsSeconds, frame)
utc = datetime(1980,1,6) + seconds(gpsSeconds - frame.gpsMinusUtcS);
if strcmp(frame.eopSource, 'iers')
    mjd = juliandate(utc) - 2400000.5;
    dut1 = deltaUT1(mjd);
    polar = polarMotion(mjd);
elseif strcmp(frame.eopSource, 'manual')
    dut1 = frame.deltaUt1S;
    polar = frame.polarMotionRad;
else
    error('LuGRE:EopSource', 'eopSource 必须是 iers 或 manual。');
end
assert(all(isfinite([dut1, polar(:)'])), 'LuGRE:EopCoverage', '地球定向参数缺失。');
% TAI-GPS 恒为 19 s，故 TAI-UTC = (GPS-UTC)+19。
rotation = dcmeci2ecef('IAU-2000/2006', datevec(utc), ...
    frame.gpsMinusUtcS + 19, dut1, polar);
end
