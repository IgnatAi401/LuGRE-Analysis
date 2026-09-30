function [rotation, rate] = icrfToItrfRotation(gpsSeconds, frame)
% ICRFTOITRFROTATION ICRF/GCRF -> ITRF 旋转矩阵及其时间导数，IAU 2000/2006；时间输入 GPST。
% GCRF 与地心 ICRF 轴对齐。使用 IERS 的 UT1、极移，不用 GMST 冒充 J2000。
% gpsSeconds 可为向量，此时返回 3x3xN（批量计算比逐个调用快得多）。
if exist('dcmeci2ecef', 'file') ~= 2
    error('LuGRE:AerospaceRequired', '坐标转换需要 Aerospace Toolbox 的 dcmeci2ecef。');
end
rotation = atTime(gpsSeconds, frame);
if nargout > 1
    h = frame.derivativeStepS;
    rate = (atTime(gpsSeconds + h, frame) - atTime(gpsSeconds - h, frame)) / (2*h);
end
end

function rotation = atTime(gpsSeconds, frame)
utc = datetime(1980,1,6) + seconds(gpsSeconds(:) - frame.gpsMinusUtcS);
if strcmp(frame.eopSource, 'iers')
    mjd = juliandate(utc) - 2400000.5;
    [dut1, polar] = iersEop(mjd);
elseif strcmp(frame.eopSource, 'manual')
    dut1 = repmat(frame.deltaUt1S, numel(utc), 1);
    polar = repmat(frame.polarMotionRad(:)', numel(utc), 1);
else
    error('LuGRE:EopSource', 'eopSource 必须是 iers 或 manual。');
end
assert(all(isfinite([dut1, polar]), 'all'), 'LuGRE:EopCoverage', '地球定向参数缺失。');
% TAI-GPS 恒为 19 s，故 TAI-UTC = (GPS-UTC)+19。
rotation = dcmeci2ecef('IAU-2000/2006', datevec(utc), ...
    repmat(frame.gpsMinusUtcS + 19, numel(utc), 1), dut1, polar);
end

function [dut1, polar] = iersEop(mjd)
% deltaUT1/polarMotion 每次调用都会重新 load aeroiersdata.mat（约 25 ms/次），
% 且只按 floor(MJD) 取当天的值（不插值）。因此按整数日缓存，结果与直接调用逐位相同。
% 更新 IERS 数据（aeroReadIERSData）后需 clear functions 以清空缓存。
persistent days dut1Cache polarCache
if isempty(days)
    days = zeros(0,1); dut1Cache = zeros(0,1); polarCache = zeros(0,2);
end
day = floor(mjd(:));
missing = setdiff(unique(day), days);
if ~isempty(missing)
    days = [days; missing];
    dut1Cache = [dut1Cache; reshape(deltaUT1(missing), [], 1)];
    polarCache = [polarCache; reshape(polarMotion(missing), [], 2)];
end
[~, index] = ismember(day, days);
dut1 = dut1Cache(index);
polar = polarCache(index,:);
end
