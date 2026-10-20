function plotErrorBars(folder)
% PLOTERRORBARS 用柱状图画 compareMethods 每个 OP 各方法相对真值的三维误差中位数（自身历元，km，对数坐标），
% 每个 OP 一张图存到 folder/figures/errorBars/<OP>.png，再把所有有结果的 OP 合成一张大图
% folder/figures/allOpErrorBars.png。方法按大类分组（参考、LS、DPE、卡尔曼滤波），同类用同一色系的深浅区分，
% 组间留空并在横轴标出类名。各图纵轴范围相同，便于跨 OP 比较；没有结果的方法在底部画灰色 ×，
% 不收敛的批处理在柱顶标 *。
root = fileparts(fileparts(mfilename('fullpath')));
if nargin < 1, folder = fullfile(root,'results','method_compare_all'); end
outFolder = fullfile(folder,'figures','errorBars');
if ~exist(outFolder,'dir'), mkdir(outFolder); end
% 大类：名称、色系（由深到浅，与组内方法一一对应）、组内方法 {字段名, 显示名}
groups = struct('name',{'参考','LS 类','DPE 类','卡尔曼滤波类'}, ...
    'colors',{{'C8475A','EBA3AC'}, {'3F4A5A','7A8699','B7BFCC'}, ...
              {'0F6B4F','2A9D74','6CC4A1','B4E2CF'}, {'0B3C5D','1F6FB2','4E93D1','8DBDE6','C6DDF2'}}, ...
    'methods',{{'nav','星上 NAV'; 'orbitPrediction','轨道预测'}, ...
               {'batchLs','批处理 LS'; 'epochLs','逐历元 LS'; 'snapshotLs','快照 LS'}, ...
               {'batchDpeGrid','批处理网格 DPE'; 'batchDpeIrls','批处理 DPE (IRLS)'; ...
                'epochDpeOrbit','逐历元网格 DPE (先验)'; 'epochDpeFree','逐历元网格 DPE (无先验)'}, ...
               {'ekfLunarStatic','EKF 月面静止'; 'ekfOrbit','EKF 动力学'; 'ekfOrbitFixedInit','EKF 动力学 (固定初值)'; ...
                'ekfCvBlend','EKF 匀速+加权'; 'ekfCvBlendFixedInit','EKF 匀速+加权 (固定初值)'}});
% 展平成方法表，并给出横轴位置（组间多空 0.8 个柱宽）
methods = {}; colors = []; x = []; groupOf = []; pos = 0;
for g = 1:numel(groups)
    m = groups(g).methods;
    methods = [methods; m]; %#ok<AGROW>
    colors = [colors; cell2mat(cellfun(@(h){sscanf(h,'%2x%2x%2x',[1 3])/255}, groups(g).colors'))]; %#ok<AGROW>
    x = [x, pos + (1:size(m,1))]; %#ok<AGROW>
    groupOf = [groupOf, repmat(g,1,size(m,1))]; %#ok<AGROW>
    pos = x(end) + 0.8;
end
groupCenters = arrayfun(@(g) mean(x(groupOf==g)), 1:numel(groups));

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
style = struct('x',x,'colors',colors,'groupCenters',groupCenters,'yLimits',yLimits);
style.groupNames = {groups.name};

% 单个 OP 的图：带数值标注，图例在右侧按大类分组。
for k = 1:numel(ops)
    figure('Visible','off','Position',[100 100 1300 560],'Color','w','Theme','light');
    layout = tiledlayout(1, 4, 'TileSpacing','compact', 'Padding','compact');
    ax = nexttile(layout, [1 3]);
    drawBars(ax, ops(k), style, true);
    ylabel(ax,'三维误差中位数 (km，对数)');
    title(ax, sprintf('%s：各方法相对真值的三维误差中位数（%s）', ...
        strrep(ops(k).window,'_','\_'), describe(ops(k).info)));
    drawLegend(nexttile(layout), groups, colors);
    exportgraphics(gcf, fullfile(outFolder,[ops(k).window '.png']), 'Resolution', 150);
    close(gcf);
end

% 合成大图：每个 OP 一格，共用图例。
nCols = 4; nRows = ceil((numel(ops)+1)/nCols); % 多留一格给图例
figure('Visible','off','Position',[50 50 1900 400*nRows],'Color','w','Theme','light');
layout = tiledlayout(nRows, nCols, 'TileSpacing','compact', 'Padding','compact');
for k = 1:numel(ops)
    ax = nexttile(layout);
    drawBars(ax, ops(k), style, false);
    ax.FontSize = 10;
    title(ax, {strrep(ops(k).window,'_','\_'), describe(ops(k).info)}, 'FontSize', 11);
    if mod(k-1,nCols) == 0, ylabel(ax,'三维误差中位数 (km)'); end
end
drawLegend(nexttile(layout), groups, colors);
title(layout, sprintf(['各 OP 各方法相对真值的三维误差中位数（自身历元，对数坐标；灰色 × 为该方法不可用或未运行）\n' ...
    '无 ≥4 星历元、所有方法均无结果的 OP 未画：%s'], strjoin(strrep(skipped,'_','\_'),'、')), 'FontSize', 14);
exportgraphics(gcf, fullfile(folder,'figures','allOpErrorBars.png'), 'Resolution', 150);
close(gcf);
fprintf('单个 OP 柱状图：%s\n合成大图：%s\n', outFolder, fullfile(folder,'figures','allOpErrorBars.png'));
end

function drawBars(ax, op, style, withLabels)
hold(ax,'on'); grid(ax,'on'); ax.YScale = 'log'; ax.YMinorGrid = 'off'; ax.XGrid = 'off';
yLimits = style.yLimits;
for j = 1:numel(style.x)
    xj = style.x(j); v = op.medianKm(j);
    h = bar(ax, xj, v, 0.85, 'FaceColor', style.colors(j,:), 'EdgeColor','none', 'BaseValue', yLimits(1));
    h.BaseLine.Visible = 'off';
    if isnan(v)
        plot(ax, xj, yLimits(1)*1.6, 'x', 'Color', [0.6 0.6 0.6], 'MarkerSize', 7, 'LineWidth', 1.2);
    elseif withLabels
        text(ax, xj, v*1.15, formatKm(v), 'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom', 'FontSize', 9);
    end
    if op.notConverged(j) && ~isnan(v)
        text(ax, xj, v*(1.15+2*withLabels), '*', 'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom', 'FontSize', 14, 'Color', 'r');
    end
end
xlim(ax, [style.x(1)-0.7, style.x(end)+0.7]); ylim(ax, yLimits);
ax.XTick = style.groupCenters; ax.XTickLabel = style.groupNames; ax.TickLength = [0 0];
end

function drawLegend(ax, groups, colors)
% 图例按大类分组：每组先放一个加粗的类名（无图标），再列组内方法。
hold(ax,'on'); axis(ax,'off');
handles = gobjects(0); labels = {}; j = 0;
for g = 1:numel(groups)
    handles(end+1) = plot(ax, NaN, NaN, 'LineStyle','none', 'Marker','none'); %#ok<AGROW>
    labels{end+1} = ['\bf' groups(g).name]; %#ok<AGROW>
    for i = 1:size(groups(g).methods,1)
        j = j + 1;
        h = bar(ax, 0, NaN, 'FaceColor', colors(j,:), 'EdgeColor','none');
        h.BaseLine.Visible = 'off';
        handles(end+1) = h; labels{end+1} = groups(g).methods{i,2}; %#ok<AGROW>
    end
end
legend(ax, handles, labels, 'Location','west', 'FontSize', 11, 'Box','off', 'Interpreter','tex');
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
