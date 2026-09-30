function files = selectRinexFiles(cfg, gpsSeconds)
% SELECTRINEXFILES 按观测 GPST 日期选择 IGS 每日广播星历，缺文件时明确报错。
% cfg.data.rinex = 'auto' 时自动匹配；否则为单个路径或路径元胞数组，原样返回。
if ~(ischar(cfg.data.rinex) || isstring(cfg.data.rinex)) || ~strcmpi(cfg.data.rinex,'auto')
    files = cellstr(cfg.data.rinex);
    return;
end
validateattributes(gpsSeconds,{'double'},{'vector','finite','nonempty'});
epoch = datetime(1980,1,6);
days = unique(dateshift(epoch + seconds([min(gpsSeconds); max(gpsSeconds)]),'start','day'));
days = (days(1):caldays(1):days(end))';
files = cell(numel(days),1);
for k = 1:numel(days)
    name = sprintf('BRDC00IGS_R_%04d%03d0000_01D_MN.rnx', year(days(k)), day(days(k),'dayofyear'));
    files{k} = fullfile(cfg.data.rinexDir, name);
    assert(isfile(files{k}), 'LuGRE:RinexMissing', ...
        '缺少 %s 的广播星历 %s（请从 IGS/CDDIS 下载并放入 %s）。', ...
        string(days(k),'yyyy-MM-dd'), name, cfg.data.rinexDir);
end
end
