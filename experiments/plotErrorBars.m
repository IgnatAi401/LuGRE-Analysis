function plotErrorBars(folder)
% PLOTERRORBARS 用柱状图画 compareMethods 每个 OP 各方法相对真值的三维误差中位数（自身历元，km，对数坐标），
% 每个 OP 一张图存到 folder/figures/errorBars/<OP>.png，再把所有有结果的 OP 合成一张大图
% folder/figures/allOpErrorBars.png。各图纵轴范围相同，便于跨 OP 比较；没有结果的方法在底部画灰色 ×，
% 不收敛的批处理在柱顶标 *。颜色沿用 plotMethodCompare。
root = fileparts(fileparts(mfilename('fullpath')));
if nargin < 1, folder = fullfile(root,'results','method_compare_all'); end
outFolder = fullfile(folder,'figures','errorBars');
if ~exist(outFolder,'dir'), mkdir(outFolder); end
methods = {'nav','星上 NAV','C8475A'; 'snapshotLs','快照 LS','C9CED6'; 'epochLs','逐历元 LS','9AA6BA'; ...
    'epochDpeOrbit','逐历元网格 DPE (先验)','8C63C9'; 'epochDpeFree','逐历元网格 DPE (无先验)','C3ADE8'; ...
    'orbitPrediction','轨道预测','D4A017'; 'batchLs','批处理 LS','5D6E8F'; ...
    'batchDpeIrls','批处理 DPE (IRLS)','2A9D8F'; 'batchDpeGrid','批处理网格 DPE','E39B2D'; ...
    'ekfOrbit','EKF 动力学','1F6FB2'; 'ekfCvBlend','EKF 匀速+加权','74A9DB'; ...
    'ekfOrbitFixedInit','EKF 动力学 (固定初值)','F28E2B'; 'ekfCvBlendFixedInit','EKF 匀速+加权 (固定初值)','F6C08A'; ...
    'ekfLunarStatic','EKF 月面静止','0B3C5D'};
colors = cell2mat(cellfun(@(h){sscanf(h,'%2x%2x%2x',[1 3])/255}, methods(:,3)));
files = dir(fullfile(folder,'OP*_*.mat'));
files = files(~contains({files.name},'_epochL1'));
key = cellfun(@(n)sscanf(n,'OP%d_%d')'*[100;1], {files.name});
[~,order] = sort(key); files = files(order);

% 读取各 OP 的中位数（km）与收敛标记；所有方法都没有结果的 OP 不画。
ops = struct('window',{},'info',{},'medianKm',{},'notConverged',{});
skipped = {};
for k = 1:numel(files)
    s = load(fullfile(folder,files(k).name),'info','stats');
    value = nan(1,size(methods,1)); flag = false(1,size(methods,1));
    for j = 1:size(methods,1)
        hit = strcmp(s.stats.method, methods{j,1});
        if any(hit), value(j) = s.stats.medianM(hit)/1000; end
        flag(j) = notConverged(s.info, methods{j,1});
    end
    if all(isnan(value)), skipped{end+1} = s.info.window; continue; end %#ok<AGROW>
    ops(end+1) = struct('window',s.info.window,'info',s.info,'medianKm',value,'notConverged',flag); %#ok<AGROW>
end
allValues = [ops.medianKm];
yLimits = 10.^[floor(log10(min(allValues))), ceil(log10(max(allValues)))+0.3];

% 单个 OP 的图：带方法名和数值标注。
for k = 1:numel(ops)
    figure('Visible','off','Position',[100 100 1100 560],'Color','w','Theme','light');
    ax = axes;
    drawBars(ax, ops(k), colors, yLimits, true);
    ax.XTick = 1:size(methods,1); ax.XTickLabel = methods(:,2); ax.XTickLabelRotation = 35;
    ylabel(ax,'三维误差中位数 (km，对数)');
    title(ax, sprintf('%s：各方法相对真值的三维误差中位数（%s）', ...
        strrep(ops(k).window,'_','\_'), describe(ops(k).info)));
    exportgraphics(gcf, fullfile(outFolder,[ops(k).window '.png']), 'Resolution', 150);
    close(gcf);
end

% 合成大图：每个 OP 一格，共用图例。
nCols = 4; nRows = ceil((numel(ops)+1)/nCols); % 多留一格给图例
figure('Visible','off','Position',[50 50 1900 400*nRows],'Color','w','Theme','light');
layout = tiledlayout(nRows, nCols, 'TileSpacing','compact', 'Padding','compact');
for k = 1:numel(ops)
    ax = nexttile(layout);
    drawBars(ax, ops(k), colors, yLimits, false);
    ax.XTick = []; ax.FontSize = 10;
    title(ax, {strrep(ops(k).window,'_','\_'), describe(ops(k).info)}, 'FontSize', 11);
    if mod(k-1,nCols) == 0, ylabel(ax,'三维误差中位数 (km)'); end
end
ax = nexttile(layout); hold(ax,'on'); axis(ax,'off');
handles = arrayfun(@(j) bar(ax, j, NaN, 'FaceColor', colors(j,:), 'EdgeColor','none'), 1:size(colors,1));
set([handles.BaseLine], 'Visible', 'off');
legend(ax, handles, methods(:,2), 'Location','west', 'FontSize', 11, 'Box','off');
title(layout, sprintf(['各 OP 各方法相对真值的三维误差中位数（自身历元，对数坐标；灰色 × 为该方法不可用或未运行）\n' ...
    '无 ≥4 星历元、所有方法均无结果的 OP 未画：%s'], strjoin(strrep(skipped,'_','\_'),'、')), 'FontSize', 14);
exportgraphics(gcf, fullfile(folder,'figures','allOpErrorBars.png'), 'Resolution', 150);
close(gcf);
fprintf('单个 OP 柱状图：%s\n合成大图：%s\n', outFolder, fullfile(folder,'figures','allOpErrorBars.png'));
end

function handles = drawBars(ax, op, colors, yLimits, withLabels)
hold(ax,'on'); grid(ax,'on'); ax.YScale = 'log'; ax.YMinorGrid = 'off';
n = size(colors,1); handles = gobjects(1,n);
for j = 1:n
    handles(j) = bar(ax, j, op.medianKm(j), 0.8, 'FaceColor', colors(j,:), 'EdgeColor','none', ...
        'BaseValue', yLimits(1));
    if isnan(op.medianKm(j))
        plot(ax, j, yLimits(1)*1.6, 'x', 'Color', [0.6 0.6 0.6], 'MarkerSize', 7, 'LineWidth', 1.2);
    elseif withLabels
        text(ax, j, op.medianKm(j)*1.15, formatKm(op.medianKm(j)), 'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom', 'FontSize', 9);
    end
    if op.notConverged(j) && ~isnan(op.medianKm(j))
        text(ax, j, op.medianKm(j)*(1.15+2*withLabels), '*', 'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom', 'FontSize', 14, 'Color', 'r');
    end
end
xlim(ax, [0.4 n+0.6]); ylim(ax, yLimits);
end

function out = describe(info)
out = sprintf('%.1f 地球半径，≥4 星历元 %d', info.meanRangeRe, info.epochsWith4);
if isfield(info,'batchDynamics') && strcmp(info.batchDynamics,'lunarStatic')
    out = [out '，批处理为月面静止模型'];
end
end

function out = formatKm(v)
if v < 10, out = sprintf('%.2f', v); elseif v < 1000, out = sprintf('%.1f', v); else, out = sprintf('%.0f', v); end
end

function flag = notConverged(info, method)
% 与 plotMethodCompare 相同：早期运行的 statusTable 对未收敛批处理也记为 ok，统一按 info.converged 判断。
map = struct('batchDpeGrid','dpeGrid','batchDpeIrls','dpeIrls','batchLs','ls');
flag = isfield(map,method) && isfield(info,'converged') && ~info.converged.(map.(method));
end
