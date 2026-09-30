function raw = readTelemetryRaw(filePath)
% READTELEMETRYRAW 读取 RAW 遥测文本：每条测量一个结构体元素。
% 基于 Simone Zocca 的 TxtParserRaw。
fid = fopen(filePath);
assert(fid > 0,'LuGRE:DataPath','无法打开 %s。',filePath);
closeFile = onCleanup(@()fclose(fid));
line = fgetl(fid);
index = 1;
while ischar(line) && ~isempty(line) % 与原解析器一致：空行即文件结束
    header = sscanf(line, "senderId: %d messageType: RAW rxTime: %f");
    parts = regexp(line, '\[ ', 'split');
    values = split(parts(2), ' ');
    values = values(2:2:end);
    count = floor(numel(values) / 7);
    for k = 1:count
        base = (k - 1) * 7;
        record = struct('rxTime',header(2),'numMeas',count, ...
            'signalId',str2double(values(base + 4)),'svId',str2double(values(base + 1)), ...
            'fdRaw',str2double(values(base + 5)),'fdRateRaw',str2double(values(base + 7)), ...
            'carrierPhase',str2double(values(base + 6)),'prRaw',str2double(values(base + 2)), ...
            'cn0',str2double(values(base + 3)));
        raw(index) = record; %#ok<AGROW> 行数未知，与原解析器一样逐条追加
        index = index + 1;
    end
    line = fgetl(fid);
end
end
