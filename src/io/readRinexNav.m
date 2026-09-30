function [gps, galileo] = readRinexNav(fileName)
% READRINEXNAV 读取 RINEX 3 广播星历中的 GPS 与 Galileo 记录。
% 每条记录保留 RINEX 字段名；TocTotalSeconds / ToeTotalSeconds 为 GPST 总秒数，
% Toe 已放到离 Toc 最近的一周（处理周日边界）。
fid = fopen(fileName, 'r');
if fid == -1
    error('LuGRE:RinexOpen', '无法打开文件: %s', fileName);
end
closeFile = onCleanup(@()fclose(fid));

blockSize = 5000; % 预分配记录数，不足时按块扩展
template = struct( ...
    'PRN', 0, 'Time', datetime, 'TocTotalSeconds', 0, 'ToeTotalSeconds', 0, ...
    'ClockBias', 0, 'ClockDrift', 0, 'ClockDriftRate', 0, ...
    'IODE', 0, 'Crs', 0, 'DeltaN', 0, 'M0', 0, ...
    'Cuc', 0, 'e', 0, 'Cus', 0, 'SqrtA', 0, ...
    'Toe', 0, 'Cic', 0, 'Omega0', 0, 'Cis', 0, ...
    'i0', 0, 'Crc', 0, 'omega', 0, 'OmegaDot', 0, ...
    'IDOT', 0, 'CodesL2', 0, 'Week', 0, 'L2PFlag', 0, ...
    'SVAccuracy', 0, 'SVHealth', 0, 'TGD', 0, 'IODC', 0, ...
    'TxTime', 0);
gps = repmat(template, blockSize, 1);
galileo = repmat(template, blockSize, 1);
gpsCount = 0;
galileoCount = 0;
gpsEpoch = datetime(1980, 1, 6);

while ~feof(fid)
    if contains(fgetl(fid), 'END OF HEADER'), break; end
end
while ~feof(fid)
    line = fgetl(fid);
    if ~ischar(line) || length(line) < 2, continue; end
    switch line(1)
        case 'G'
            gpsCount = gpsCount + 1;
            if gpsCount > numel(gps), gps = [gps; repmat(template, blockSize, 1)]; end %#ok<AGROW>
            gps(gpsCount) = parseRecord(readBlock(fid, line, 7), gpsEpoch);
        case 'E'
            galileoCount = galileoCount + 1;
            if galileoCount > numel(galileo), galileo = [galileo; repmat(template, blockSize, 1)]; end %#ok<AGROW>
            galileo(galileoCount) = parseRecord(readBlock(fid, line, 7), gpsEpoch);
        case {'C', 'J', 'I'}
            readBlock(fid, line, 7);
        case {'R', 'S'}
            readBlock(fid, line, 3);
    end
end
gps = gps(1:gpsCount);
galileo = galileo(1:galileoCount);
fprintf('读取 %s：GPS %d 条，Galileo %d 条。\n', fileName, gpsCount, galileoCount);
end

function block = readBlock(fid, firstLine, dataLineCount)
block = cell(1 + dataLineCount, 1);
block{1} = firstLine;
for k = 1:dataLineCount
    block{1+k} = fgetl(fid);
end
end

function record = parseRecord(lines, gpsEpoch)
first = lines{1};
record.PRN = str2double(first(2:3));
year = str2double(first(5:8));
month = str2double(first(10:11));
day = str2double(first(13:14));
hour = str2double(first(16:17));
minute = str2double(first(19:20));
second = str2double(first(22:23));
record.Time = datetime(year, month, day, hour, minute, second);
record.TocTotalSeconds = round(seconds(record.Time - gpsEpoch));
record.ToeTotalSeconds = NaN;
record.ClockBias = parseNumber(first(24:42));
record.ClockDrift = parseNumber(first(43:61));
record.ClockDriftRate = parseNumber(first(62:80));

line = pad(lines{2}, 80);
record.IODE = parseNumber(line(5:23));
record.Crs = parseNumber(line(24:42));
record.DeltaN = parseNumber(line(43:61));
record.M0 = parseNumber(line(62:80));

line = pad(lines{3}, 80);
record.Cuc = parseNumber(line(5:23));
record.e = parseNumber(line(24:42));
record.Cus = parseNumber(line(43:61));
record.SqrtA = parseNumber(line(62:80));

line = pad(lines{4}, 80);
record.Toe = parseNumber(line(5:23));
record.ToeTotalSeconds = record.Toe + 604800 * round((record.TocTotalSeconds - record.Toe)/604800);
record.Cic = parseNumber(line(24:42));
record.Omega0 = parseNumber(line(43:61));
record.Cis = parseNumber(line(62:80));

line = pad(lines{5}, 80);
record.i0 = parseNumber(line(5:23));
record.Crc = parseNumber(line(24:42));
record.omega = parseNumber(line(43:61));
record.OmegaDot = parseNumber(line(62:80));

line = pad(lines{6}, 80);
record.IDOT = parseNumber(line(5:23));
record.CodesL2 = parseNumber(line(24:42));
record.Week = parseNumber(line(43:61));
record.L2PFlag = parseNumber(line(62:80));

line = pad(lines{7}, 80);
record.SVAccuracy = parseNumber(line(5:23));
record.SVHealth = parseNumber(line(24:42));
record.TGD = parseNumber(line(43:61));
record.IODC = parseNumber(line(62:80));

line = pad(lines{8}, 80);
record.TxTime = parseNumber(line(5:23));
end

function value = parseNumber(text)
% RINEX 数值可能用 D 作指数；空字段或无法解析时记为 0。
if isempty(text)
    value = 0; return;
end
value = str2double(strrep(text, 'D', 'E'));
if isnan(value), value = 0; end
end
