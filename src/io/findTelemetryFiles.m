function [rawPath, navPath] = findTelemetryFiles(dataRoot, window)
% FINDTELEMETRYFILES 定位观测窗口的 RAW/NAV 遥测文本（L0/TLM/TLM_{RAW,NAV}_*_<window>.txt）。
% RAW 必须唯一存在；NAV 缺失时返回 ''，由调用者决定是否必需。
folder = fullfile(dataRoot,'L0','TLM');
raw = dir(fullfile(folder,['TLM_RAW_*_' window '.txt']));
assert(isscalar(raw),'LuGRE:DataPath','窗口 %s 的 RAW 文件数为 %d（应为 1）。',window,numel(raw));
rawPath = fullfile(folder,raw.name);
nav = dir(fullfile(folder,['TLM_NAV_*_' window '.txt']));
navPath = '';
if isscalar(nav), navPath = fullfile(folder,nav.name); end
end
