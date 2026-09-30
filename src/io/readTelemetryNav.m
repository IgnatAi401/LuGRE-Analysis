function nav = readTelemetryNav(filePath)
% READTELEMETRYNAV 读取星上 NAV 解遥测文本：每个历元一个结构体元素。
% 基于 Simone Zocca 的 TxtParserNav。
fid = fopen(filePath);
assert(fid > 0,'LuGRE:DataPath','无法打开 %s。',filePath);
closeFile = onCleanup(@()fclose(fid));
format = "senderId: %d messageType: NAV " + ...
    "rxTime: %f appName: NAV wn: %d tow: %d decimals: %f nSat: %d " + ...
    "posX: %f posY: %f posZ: %f " + ...
    "velX: %f velY: %f velZ: %f " + ...
    "posStd: %f velStd: %f timStd: %f " + ...
    "clockBias: %f clockDrift: %f ggto: %f " + ...
    "GDOP: %f PDOP: %f HDOP: %f VDOP: %f TDOP: %f";
names = {'posX','posY','posZ','velX','velY','velZ','posStd','velStd','timStd', ...
    'clockBias','clockDrift','ggto','GDOP','PDOP','HDOP','VDOP','TDOP'};
line = fgetl(fid);
index = 1;
while ischar(line) && ~isempty(line) % 与原解析器一致：空行即文件结束
    values = sscanf(line, format);
    record = struct('rxTime',values(2),'wn',values(3),'tow',values(4) + values(5),'nSat',values(6));
    for k = 1:numel(names)
        record.(names{k}) = values(6 + k);
    end
    nav(index) = record; %#ok<AGROW> 行数未知，逐条追加
    line = fgetl(fid);
    index = index + 1;
end
end
