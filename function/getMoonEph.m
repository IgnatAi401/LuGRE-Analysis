function results = getMoonEph(gpsSeconds, cfg)
% GETMOONEPH 地心几何月球星历：ICRF 赤道轴，UTC 时间标签。
% 显式传 cfg：返回每个输入时刻的 [JD_UTC, X,Y,Z (km), VX,VY,VZ (km/s)]。
% 仅传时间：兼容历史脚本，返回从起点开始的分钟采样历表。
% 保留旧接口单位；只应将第 2:7 列乘 1000，禁止转换 JD 列。
legacySampling = nargin < 2;
if legacySampling
    cfg = main5_config();
end
validateattributes(gpsSeconds, {'double'}, {'vector','finite','nonempty'});
gpsSeconds = gpsSeconds(:);
epoch = datetime(1980,1,6);
% 前后各扩展一个采样间隔，覆盖小数秒及单历元请求。
step = cfg.moon.stepMinutes * 60;
startGps = floor(min(gpsSeconds)) - step;
stopGps = ceil(max(gpsSeconds)) + step;
startUtc = epoch + seconds(startGps - cfg.frame.gpsMinusUtcS);
stopUtc = epoch + seconds(stopGps - cfg.frame.gpsMinusUtcS);
response = webread(cfg.moon.url, ...
    'format','json', 'COMMAND','''301''', 'CENTER','''500@399''', ...
    'OBJ_DATA','''NO''', 'MAKE_EPHEM','''YES''', 'EPHEM_TYPE','''VECTORS''', ...
    'REF_SYSTEM','''ICRF''', 'REF_PLANE','''FRAME''', 'TIME_TYPE','''UT''', ...
    'VEC_CORR','''NONE''', 'OUT_UNITS','''KM-S''', 'VEC_TABLE','''2''', ...
    'CSV_FORMAT','''YES''', ...
    'START_TIME',['''' char(string(startUtc,'yyyy-MM-dd HH:mm:ss')) ''''], ...
    'STOP_TIME',['''' char(string(stopUtc,'yyyy-MM-dd HH:mm:ss')) ''''], ...
    'STEP_SIZE',sprintf('''%dm''', cfg.moon.stepMinutes), ...
    weboptions('Timeout',cfg.moon.timeoutS));
if ~isfield(response,'result')
    error('LuGRE:HorizonsResponse', 'Horizons 未返回星历：%s', jsonencode(response));
end
raw = response.result;
first = strfind(raw,'$$SOE'); last = strfind(raw,'$$EOE');
if numel(first) ~= 1 || numel(last) ~= 1
    error('LuGRE:HorizonsResponse', 'Horizons 星历段缺失：%s', raw);
end
lines = splitlines(strtrim(string(raw(first+5:last-1))));
table = nan(numel(lines),7);
for k = 1:numel(lines)
    values = split(strtrim(lines(k)), ',');
    if numel(values) < 8
        error('LuGRE:HorizonsParse', '无法解析月球星历第 %d 行。', k);
    end
    table(k,:) = str2double(values([1,3:8]))';
end
assert(all(isfinite(table),'all') && all(diff(table(:,1)) > 0), ...
    'LuGRE:HorizonsParse', '月球星历包含无效数值或时间未递增。');
% 在小量相对时间轴插值，保留真实时间间隔；不再按 NAV 行号匹配。
sampleSeconds = (table(:,1) - table(1,1))*86400;
firstGps = (table(1,1) - 2444244.5)*86400 + cfg.frame.gpsMinusUtcS;
if legacySampling
    sampleGps = firstGps + sampleSeconds;
    keep = sampleGps >= floor(min(gpsSeconds))-1e-3 & ...
        sampleGps <= ceil(max(gpsSeconds))+1e-3;
    results = table(keep,:);
    return;
end
state = interp1(sampleSeconds, table(:,2:7), gpsSeconds-firstGps, 'pchip');
assert(all(isfinite(state),'all'), 'LuGRE:MoonCoverage', '月球星历未覆盖观测时间。');
results = [2444244.5 + (gpsSeconds-cfg.frame.gpsMinusUtcS)/86400, state];
end
