function plotLunarTopView(folder)
% PLOTLUNARTOPVIEW 月面静止窗口的定位结果俯视图：以真值着陆点为原点，在当地东-北水平面上画各方法
% 定位点相对着陆点的位置，颜色表示高度（当地天向）误差。每个月面 OP 一张图存到
% results/lunar_static/topView/<OP>.png，另存所有 OP 的总览 allOps.png 和数字汇总 topViewSummary.csv。
% 当地东-北-天坐标在月固系（DE405 主轴系）中定义于真值中位点；误差向量 e = 估计 - 真值（ECEF）
% 逐历元旋转到月固系 R_mf·R_itrf'·e（同一时刻的差向量只需旋转，不用月心星历）。
% 批处理解直接用月固系 p 减真值。图中灰色虚线为地球方向在水平面的投影（误差主要沿“月→地”视线）。
% 背景为着陆点周围的 LRO WAC 月面影像（src/plot/lunarBasemap.m，瓦片缓存在 results/lunar_basemap/）。
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'src')));
if nargin < 1, folder = fullfile(root,'results','method_compare_all'); end
outFolder = fullfile(root,'results','lunar_static','topView');
if ~exist(outFolder,'dir'), mkdir(outFolder); end
tileFolder = fullfile(root,'results','lunar_basemap');
cfg = lugreConfig();
truthAll = load(fullfile(root,'results','lunar_static','truthMoonFixed.mat'));
windows = {'OP38_0','OP40_0','OP74_0','OP76_0','OP77_0','OP77_1','OP78_1'};
clouds = {'epochLs','逐历元 LS'; 'nav','星上 NAV'; 'ekfLunarStatic','EKF 月面静止'};
batches = {'batchDpeGrid','dpeGrid','批处理网格 DPE','p'; 'batchDpeIrls','dpeIrls','批处理 DPE (IRLS)','d'; ...
    'batchLs','ls','批处理 LS','s'};
cmap = divergingMap();

data = cell(numel(windows),1);
for w = 1:numel(windows)
    s = load(fullfile(folder,[windows{w} '.mat']),'series','info');
    truthMf = truthAll.(windows{w}).median(:);
    [east, north, up] = enuBasis(truthMf);
    basis = [east north up]'; % 月固系 -> ENU
    % 需要旋转的历元：任一点云方法有值且真值有效
    used = false(numel(s.series.gpsSeconds),1);
    for j = 1:size(clouds,1)
        if isfield(s.series,clouds{j,1}), used = used | isfinite(s.series.(clouds{j,1}).distanceM); end
    end
    used = used & all(isfinite(s.series.truth),2);
    t = s.series.gpsSeconds(used);
    toMoonFixed = pagemtimes(icrfToMoonFixedRotation(t), 'none', icrfToItrfRotation(t, cfg.frame), 'transpose');
    rotate = @(v) permute(pagemtimes(toMoonFixed, permute(v,[2 3 1])),[3 1 2]) * basis'; % n x 3 ECEF -> ENU
    d.window = windows{w};
    d.truthMf = truthMf;
    earth = rotate(-s.series.truth(used,:)./vecnorm(s.series.truth(used,:),2,2));
    d.earthAzimuth = atan2d(mean(earth(:,1)), mean(earth(:,2))); % 从北向东
    d.earthElevation = asind(mean(earth(:,3)));
    for j = 1:size(clouds,1)
        name = clouds{j,1};
        if ~isfield(s.series,name), d.(name) = zeros(0,3); continue; end
        e = rotate(s.series.(name).positionEcef(used,:) - s.series.truth(used,:));
        d.(name) = e(all(isfinite(e),2),:)/1000;
    end
    for j = 1:size(batches,1)
        d.(batches{j,1}) = (basis*(s.info.moonFixedState.(batches{j,2}) - truthMf))'/1000;
    end
    data{w} = d;
    drawWindow(d, clouds, batches, cmap, fullfile(outFolder,[windows{w} '.png']), tileFolder);
end
drawOverview(data, batches, cmap, fullfile(outFolder,'allOps.png'), tileFolder);
% 汇总（km、度）：地球方位/仰角，各点云历元数与水平误差中位数、EKF 末历元，三种批处理 ENU 误差
summary = table();
for w = 1:numel(data)
    d = data{w};
    row = table(string(d.window), d.earthAzimuth, d.earthElevation, 'VariableNames', {'window','earthAzimuthDeg','earthElevationDeg'});
    for j = 1:size(clouds,1)
        p = d.(clouds{j,1});
        row.([clouds{j,1} 'Count']) = size(p,1);
        row.([clouds{j,1} 'HorizontalMedianKm']) = median(vecnorm(p(:,1:2),2,2));
    end
    last = nan(1,3); if ~isempty(d.ekfLunarStatic), last = d.ekfLunarStatic(end,:); end
    row = [row array2table(last, 'VariableNames', {'ekfLastEastKm','ekfLastNorthKm','ekfLastUpKm'})]; %#ok<AGROW>
    for j = 1:size(batches,1)
        row = [row array2table(d.(batches{j,1}), 'VariableNames', strcat(batches{j,2}, {'EastKm','NorthKm','UpKm'}))]; %#ok<AGROW>
    end
    summary = [summary; row]; %#ok<AGROW>
end
writetable(summary, fullfile(outFolder,'topViewSummary.csv'));
disp(summary(:,{'window','earthAzimuthDeg','earthElevationDeg','dpeGridEastKm','dpeGridNorthKm','dpeGridUpKm'}));
fprintf('图已保存到 %s\n', outFolder);
end

function drawWindow(d, clouds, batches, cmap, file, tileFolder)
% 单个 OP（2×2）：逐历元 LS、NAV、EKF 月面静止三种逐历元点云（前两种叠加批处理解），
% 第四栏放大到三种批处理解（每种一个固定点）。
figure('Visible','off','Position',[100 100 1500 1250],'Color','w','Theme','light');
layout = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
title(layout, sprintf('%s 月面静止：定位点相对着陆点的俯视位置（颜色 = 高度误差）', ...
    strrep(d.window,'_','\_')), 'FontSize', 18, 'FontWeight', 'bold');
for j = 1:4
    ax = nexttile(layout); hold(ax,'on'); grid(ax,'on'); box(ax,'on'); ax.FontSize = 15;
    if j <= 3
        points = d.(clouds{j,1});
        if isempty(points)
            axis(ax,'off'); title(ax,[clouds{j,2} '：无结果（没有 ≥5 星历元，无法初始化）'],'FontWeight','normal');
            continue;
        end
        limit = robustLimit(points(:,1:2), 0.98);
        outside = nnz(max(abs(points(:,1:2)),[],2) > limit);
        drawCloud(ax, points, 6);
        subtitle = sprintf('%s（%d 个历元，水平误差中位数 %.1f km）', clouds{j,2}, size(points,1), ...
            median(vecnorm(points(:,1:2),2,2)));
        if j == 3, subtitle = sprintf('%s\n末历元：东 %+.1f、北 %+.1f、高 %+.1f km', subtitle, points(end,:)); end
        if outside > 0, subtitle = sprintf('%s\n另有 %d 点（%.1f%%）在框外', subtitle, outside, ...
                100*outside/size(points,1)); end
        cl = robustLimit(points(:,3), 0.95);
    else
        batchPoints = cell2mat(cellfun(@(f) d.(f), batches(:,1), 'UniformOutput', false));
        limit = 1.3*max(abs(batchPoints(:,1:2)),[],'all');
        cl = max(abs(batchPoints(:,3)));
        subtitle = '批处理解（月固系固定点，每种方法一个点）';
    end
    limit = max(limit, 0.05);
    axis(ax,'equal'); xlim(ax,[-limit limit]); ylim(ax,[-limit limit]);
    drawBasemap(ax, d, limit, tileFolder);
    drawEarth(ax, d, limit);
    plot(ax, 0, 0, '+', 'MarkerSize', 18, 'Color', 'k', 'LineWidth', 2, 'Tag', 'legend', 'DisplayName', '着陆点（真值）');
    if j ~= 3, drawBatch(ax, d, batches, 'k', j == 4); end
    colormap(ax, cmap); clim(ax, [-1 1]*max(cl, 1e-3));
    cb = colorbar(ax); cb.Label.String = '高度误差（km，正 = 偏高）';
    xlabel(ax,'东向误差 (km)'); ylabel(ax,'北向误差 (km)'); title(ax, subtitle, 'FontWeight','normal');
    if j == 4, legend(ax, flipud(findobj(ax,'Tag','legend')), 'Location','northeast', 'FontSize', 14); end
end
exportgraphics(gcf, file, 'Resolution', 150);
close(gcf);
end

function drawOverview(data, batches, cmap, file, tileFolder)
% 所有 OP：逐历元 LS 点云 + 三种批处理解，各 OP 各自的坐标范围与色标。
figure('Visible','off','Position',[100 100 1900 860],'Color','w','Theme','light');
layout = tiledlayout(2,4,'TileSpacing','loose','Padding','loose');
title(layout, '月面静止窗口：逐历元 LS（点）与批处理解（黑框标记）相对着陆点的俯视位置，颜色 = 高度误差', ...
    'FontSize', 17, 'FontWeight', 'bold');
for w = 1:numel(data)
    d = data{w};
    ax = nexttile(layout); hold(ax,'on'); grid(ax,'on'); box(ax,'on'); ax.FontSize = 13;
    points = d.epochLs;
    limit = robustLimit(points(:,1:2), 0.98);
    drawCloud(ax, points, 4);
    axis(ax,'equal'); xlim(ax,[-limit limit]); ylim(ax,[-limit limit]);
    drawBasemap(ax, d, limit, tileFolder);
    drawEarth(ax, d, limit);
    plot(ax, 0, 0, '+', 'MarkerSize', 16, 'Color', 'k', 'LineWidth', 2);
    drawBatch(ax, d, batches, 'k');
    colormap(ax, cmap); clim(ax, [-1 1]*robustLimit(points(:,3), 0.95));
    cb = colorbar(ax); cb.Label.String = '高度误差 (km)';
    xlabel(ax,'东 (km)'); ylabel(ax,'北 (km)');
    title(ax, sprintf('%s  网格 DPE：水平 %.2f、高 %+.2f km', strrep(d.window,'_','\_'), ...
        norm(d.batchDpeGrid(1:2)), d.batchDpeGrid(3)), 'FontWeight','normal');
end
ax = nexttile(layout); axis(ax,'off'); hold(ax,'on');
for j = 1:size(batches,1)
    plot(ax, nan, nan, batches{j,4}, 'MarkerSize', 10, 'MarkerFaceColor', [0.85 0.85 0.85], ...
        'MarkerEdgeColor', 'k', 'DisplayName', batches{j,3});
end
plot(ax, nan, nan, '+', 'MarkerSize', 16, 'Color', 'k', 'LineWidth', 2, 'DisplayName', '着陆点（真值）');
plot(ax, nan, nan, '--', 'Color', [0.15 0.15 0.15], 'LineWidth', 1.4, 'DisplayName', '地球方向（水平投影）');
legend(ax, 'Location', 'west', 'FontSize', 12);
text(ax, 0.02, 0.12, {'每个 OP 色标不同；点云坐标框取', '98% 分位，个别离群点在框外', '背景：LRO WAC 月面影像（NASA）'}, ...
    'Units','normalized','FontSize',11,'Color',[0.3 0.3 0.3]);
exportgraphics(gcf, file, 'Resolution', 150);
close(gcf);
end

function drawCloud(ax, points, markerSize)
% 按高度误差绝对值从小到大绘制，误差大的点在上层。
[~, order] = sort(abs(points(:,3)));
scatter(ax, points(order,1), points(order,2), markerSize, points(order,3), 'filled', 'MarkerFaceAlpha', 0.8);
end

function drawBatch(ax, d, batches, edge, withLegend)
if nargin < 5, withLegend = false; end
for j = 1:size(batches,1)
    p = d.(batches{j,1});
    h = scatter(ax, p(1), p(2), 140, p(3), 'filled', batches{j,4}, 'MarkerEdgeColor', edge, 'LineWidth', 1.3);
    if withLegend
        h.Tag = 'legend'; h.DisplayName = sprintf('%s（高度 %+.2f km）', batches{j,3}, p(3));
    end
end
end

function drawBasemap(ax, d, limit, tileFolder)
% 月面影像垫底（正射投影，月盘外留白）：拉伸到 2%~98% 分位后向白色混合，保证彩色点和黑色标记清晰；网格线画在影像上层。
[img, x, y] = lunarBasemap(d.truthMf, limit, 900, tileFolder);
outside = isnan(img);
low = prctile(img(:), 2); high = prctile(img(:), 98);
img = min(max((img - low)/max(high - low, eps), 0), 1);
img = 0.45 + 0.55*img;
img(outside) = 1; % 月盘外留白（min/max 会忽略 NaN，需先记下）
h = image(ax, x, y, repmat(img, 1, 1, 3));
uistack(h, 'bottom');
ax.YDir = 'normal'; ax.Layer = 'top'; ax.GridColor = [1 1 1]; ax.GridAlpha = 0.35;
end

function drawEarth(ax, d, limit)
direction = [sind(d.earthAzimuth) cosd(d.earthAzimuth)];
plot(ax, [0 direction(1)]*1.4*limit, [0 direction(2)]*1.4*limit, '--', 'Color', [0.15 0.15 0.15], 'LineWidth', 1.4);
text(ax, direction(1)*0.8*limit, direction(2)*0.8*limit, sprintf('地球方向\n仰角 %.0f°', d.earthElevation), ...
    'Color', [0.1 0.1 0.1], 'FontSize', 12, 'VerticalAlignment', 'bottom', 'BackgroundColor', [1 1 1], 'Margin', 1);
end

function limit = robustLimit(values, level)
limit = prctile(max(abs(values),[],2), 100*level);
end

function [east, north, up] = enuBasis(p)
lat = asin(p(3)/norm(p)); lon = atan2(p(2), p(1));
east = [-sin(lon); cos(lon); 0];
north = [-sin(lat)*cos(lon); -sin(lat)*sin(lon); cos(lat)];
up = [cos(lat)*cos(lon); cos(lat)*sin(lon); sin(lat)];
end

function cmap = divergingMap()
% 蓝（偏低）- 浅灰 - 红（偏高），颜色越深误差越大；中心用浅灰以便在白底上可见。
anchors = [33 102 172; 103 169 207; 217 217 217; 239 138 98; 178 24 43]/255;
cmap = interp1(linspace(0,1,5), anchors, linspace(0,1,256));
end
