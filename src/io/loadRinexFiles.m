function [gps, galileo] = loadRinexFiles(files)
% LOADRINEXFILES 读取并合并多个 RINEX 导航文件（跨日观测）。
files = cellstr(files);
gps = []; galileo = [];
for k = 1:numel(files)
    [g, e] = readRinexNav(files{k});
    gps = [gps; g(:)]; %#ok<AGROW>
    galileo = [galileo; e(:)]; %#ok<AGROW>
end
end
