function results = replotLugreTruth(matFile, outputDir)
% REPLOTLUGRETRUTH 无需重新定位，读取 main5 保存的 results 并保存真值图。
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root,'function'));
if nargin<1, matFile = fullfile(root,'try1.mat'); end
if nargin<2, outputDir = fullfile(root,'results','try1_truth'); end
saved = load(matFile,'results');
assert(isfield(saved,'results'),'LuGRE:SavedResults','MAT 文件缺少 results。');
cfg = main5_config();
results = compareLugreTruth(saved.results,cfg.truth);
figures = plotLugreTruth(results);
if ~isfolder(outputDir), mkdir(outputDir); end
names = {'truth_3d_error','truth_xyz_error','truth_trajectory'};
for k = 1:numel(figures)
    savefig(figures(k),fullfile(outputDir,[names{k},'.fig']));
    exportgraphics(figures(k),fullfile(outputDir,[names{k},'.png']),'Resolution',150);
end
save(fullfile(outputDir,'truth_comparison.mat'),'results');
end
