function plotMethodCompare(folder)
% PLOTMETHODCOMPARE 为 compareMethods 的每个 OP 结果画各方法三维误差随时间的曲线（对数坐标），
% 保存为 folder/figures/<OP>.png。颜色与汇报 PPT 一致；没有结果的方法不画，不收敛的批处理在图例中标注。
root = fileparts(fileparts(mfilename('fullpath')));
if nargin < 1, folder = fullfile(root,'results','method_compare_all'); end
outFolder = fullfile(folder,'figures');
if ~exist(outFolder,'dir'), mkdir(outFolder); end
methods = {'batchDpeGrid','批处理网格 DPE','E39B2D'; 'batchDpeIrls','批处理 DPE (IRLS)','2A9D8F'; ...
    'batchLs','批处理 LS','5D6E8F'; 'epochDpeOrbit','逐历元网格 DPE (先验)','8C63C9'; ...
    'epochLs','逐历元 LS','9AA6BA'; 'ekfOrbit','EKF 动力学','1F6FB2'; 'ekfCvBlend','EKF 匀速+加权','74A9DB'; ...
    'ekfLunarStatic','EKF 月面静止','0B3C5D'; 'nav','星上 NAV','C8475A'};
files = dir(fullfile(folder,'OP*_*.mat'));
for k = 1:numel(files)
    s = load(fullfile(folder,files(k).name),'series','info','statusTable');
    present = methods(isfield(s.series, methods(:,1)),:);
    if isempty(present), continue; end
    t = (s.series.gpsSeconds - s.series.gpsSeconds(1))/60; unit = 'min';
    if t(end) > 180, t = t/60; unit = 'h'; end % 月面长窗口用小时
    figure('Visible','off','Position',[100 100 1100 520],'Color','w','Theme','light');
    ax = axes; hold(ax,'on'); grid(ax,'on'); ax.YScale = 'log';
    names = {};
    for j = 1:size(present,1)
        d = s.series.(present{j,1}).distanceM/1000;
        color = sscanf(present{j,3},'%2x%2x%2x',[1 3])/255;
        if strcmp(present{j,1},'nav') || startsWith(present{j,1},'epoch')
            plot(ax, t, d, '.', 'Color', color, 'MarkerSize', 4); % 逐历元/NAV 只在 NAV 历元有值，可能有空档
        else
            plot(ax, t, d, '-', 'Color', color, 'LineWidth', 1.2);
        end
        static = startsWith(present{j,1},'batch') && isfield(s.info,'batchDynamics') && ...
            strcmp(s.info.batchDynamics,'lunarStatic');
        names{end+1} = [present{j,2} ternary(static,'（月面静止）','') ...
            ternary(notConverged(s.info,present{j,1}),' *（不收敛）','')]; %#ok<AGROW>
    end
    xlabel(ax,['窗口内时间 (' unit ')']); ylabel(ax,'三维误差 (km，对数)');
    title(ax, sprintf('%s：各方法相对真值的三维误差（距离 %.1f 地球半径，≥4 星历元 %d）', ...
        strrep(s.info.window,'_','\_'), s.info.meanRangeRe, s.info.epochsWith4));
    legend(ax, names, 'Location','eastoutside');
    exportgraphics(gcf, fullfile(outFolder,[s.info.window '.png']), 'Resolution', 150);
    close(gcf);
end
fprintf('图已保存到 %s\n', outFolder);
end

function flag = notConverged(info, method)
% 早期运行的 statusTable 对未收敛批处理也记为 ok，统一按 info.converged 判断。
map = struct('batchDpeGrid','dpeGrid','batchDpeIrls','dpeIrls','batchLs','ls');
flag = isfield(map,method) && isfield(info,'converged') && ~info.converged.(map.(method));
end

function out = ternary(condition, a, b)
if condition, out = a; else, out = b; end
end
