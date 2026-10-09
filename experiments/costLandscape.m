function land = costLandscape(window, outputFolder, replotOnly)
% COSTLANDSCAPE 真实数据的 L1 代价曲面，供汇报“远距离代价曲面是一条又长又窄的斜谷”。
% 左：逐历元。选一个代表历元（网格解距 L1 最优为中位数），在“径向偏移 - 钟差偏移”平面上
%     （横向固定在 L1 最优）画单历元 L1 代价，叠加逐历元网格的格点（钟差轴是真实网格轴）。
% 右：批处理 dpeGrid。在最终解处线性化，取白化坐标中最弱的两个方向 u1、u2（其余 4 维为 0），
%     每历元钟差取中位数（与 solveBatchGrid 的代价完全相同），叠加一级 5x5 网格。
% 颜色为代价增量 ΔJ = J - min J（m，对数色标）。数据与图存到 outputFolder。
% replotOnly=true 时只读取已存的 .mat 重新作图。
% 例：costLandscape('OP5_0')
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, genpath(fullfile(root,'src')));
if nargin < 2 || isempty(outputFolder), outputFolder = fullfile(root,'results','cost_landscape'); end
if ~exist(outputFolder,'dir'), mkdir(outputFolder); end
dataFile = fullfile(outputFolder,[window '_landscape.mat']);
if nargin >= 3 && replotOnly
    load(dataFile,'land'); plotLandscape(land, outputFolder); return;
end
cfg = lugreConfig(); cfg.data.window = window; cfg.run.plot = false; cfg.run.progressEvery = 0;
cfg.enable.orbit = false; % 与 epochL1Check 一致：只看观测代价
if strcmp(window,'OP14_0'), cfg.run.skipEpochs = 3; end
land = struct('window',window);
land.epoch = epochLandscape(cfg);
land.batch = batchLandscape(cfg);
save(dataFile,'land');
plotLandscape(land, outputFolder);
end

%% 逐历元
function e = epochLandscape(cfg)
[rawPath,navPath] = findTelemetryFiles(cfg.data.root,cfg.data.window);
nav = readTelemetryNav(navPath); nav = nav(cfg.run.skipEpochs+1:end);
raw = readTelemetryRaw(rawPath);
time = [nav.rxTime]'; n = numel(time);
[gps,galileo] = loadRinexFiles(selectRinexFiles(cfg,time));
raw = addTransmitTime(raw(ismember([raw.rxTime],time)),cfg.constants.c);
raw = addSatelliteStates(raw,gps,galileo,cfg.ephemeris.maxAgeS);
truth = loadTruthEcef(time, cfg.truth, cfg.frame);

% 与 epochL1Check 相同的链式求解，记下每历元网格中心（上一有效解）以复原格点位置。
dpePosition = [nav(1).posX; nav(1).posY; nav(1).posZ]; dpeClock = cfg.initial.clockBiasM;
lsPosition = dpePosition; lsClock = dpeClock;
rec = repmat(struct('obs',[],'center',[],'centerClock',NaN,'grid',[],'l1',[],'ls',[]), n, 1);
distance = nan(n,1);
for i = 1:n
    obs = buildObservations(selectEpochMeasurements(time(i),raw),dpePosition,cfg);
    rec(i).center = dpePosition; rec(i).centerClock = dpeClock;
    grid = solveEpochDpe(obs,dpePosition,dpeClock,nan(3,1),i,cfg);
    ls = solveEpochLs(obs,lsPosition,lsClock,cfg);
    if grid.valid, dpePosition = grid.position; dpeClock = grid.clockM; end
    if ls.valid, lsPosition = ls.position; lsClock = ls.clockM; end
    if ~(grid.valid && ls.valid && all(isfinite(truth(i,:)))), continue; end
    l1 = solveL1(obs, ls.position, ls.clockM);
    rec(i).obs = obs; rec(i).grid = grid; rec(i).l1 = l1; rec(i).ls = ls;
    distance(i) = norm(grid.position - l1.position);
end
good = find(isfinite(distance));
[~,pick] = min(abs(distance(good) - median(distance(good))));
i = good(pick); r = rec(i);

% 平面：径向单位向量 u（L1 解方向）与钟差；横向固定在 L1 最优。
p0 = r.l1.position; c0 = r.l1.clockM; u = p0/norm(p0);
cost = @(p,c) sum(abs(r.obs.range - vecnorm(p'-r.obs.satPosition,2,2) - c));
gridRadial = u'*(r.grid.position - p0); gridClock = r.grid.clockM - c0;
gridTransverse = norm((r.grid.position - p0) - gridRadial*u);
half = max(5e3, 1.3*max(abs([gridRadial gridClock])));
axisR = linspace(-half, half, 601); axisC = linspace(-half, half, 601);
[R,C] = meshgrid(axisR, axisC);
J = zeros(size(R));
for j = 1:numel(r.obs.range)
    rho = vecnorm(p0' + R(:).*u' - r.obs.satPosition(j,:),2,2);
    J(:) = J(:) + abs(r.obs.range(j) - rho - (c0 + C(:)));
end
% 窄谷放大：以网格解为中心、宽 3 个格距，看格点与谷底的关系。
step = cfg.dpe.coarseOffsetsM; if i >= cfg.dpe.switchEpoch, step = cfg.dpe.fineOffsetsM; end
spacing = step(2) - step(1);
zoomHalf = 3*spacing;
zr = gridRadial + linspace(-zoomHalf, zoomHalf, 401); zc = gridClock + linspace(-zoomHalf, zoomHalf, 401);
[ZR,ZC] = meshgrid(zr, zc); ZJ = zeros(size(ZR));
for j = 1:numel(r.obs.range)
    rho = vecnorm(p0' + ZR(:).*u' - r.obs.satPosition(j,:),2,2);
    ZJ(:) = ZJ(:) + abs(r.obs.range(j) - rho - (c0 + ZC(:)));
end
% 格点（粗网格）：钟差轴是真实网格轴；位置格点沿径向的投影间距取同一步长。
latticeR = u'*(r.center - p0) + step; latticeC = r.centerClock - c0 + cfg.dpe.clockOffsetsM;
e = struct('epochIndex',i,'gpsSeconds',time(i),'satCount',numel(r.obs.range), ...
    'radialM',axisR,'clockM',axisC,'cost',J,'costMin',cost(p0,c0), ...
    'zoomRadialM',zr,'zoomClockM',zc,'zoomCost',ZJ,'spacingM',spacing, ...
    'latticeRadialM',latticeR,'latticeClockM',latticeC, ...
    'gridRadialM',gridRadial,'gridClockM',gridClock,'gridTransverseM',gridTransverse, ...
    'gridCost',cost(r.grid.position,r.grid.clockM), ...
    'lsRadialM',u'*(r.ls.position-p0),'lsClockM',r.ls.clockM-c0, ...
    'truthRadialM',u'*(truth(i,:)'-p0),'gridToL1MedianM',median(distance(good)), ...
    'gridToL1M',distance(i),'validEpochs',numel(good));
fprintf('%s 逐历元：代表历元 %d（%d 星），网格解距 L1 最优 %.1f km（中位数 %.1f km），径向 %.1f km、钟差 %.1f km、横向 %.2f km；代价 网格 %.1f / L1 %.1f m\n', ...
    cfg.data.window, i, e.satCount, distance(i)/1e3, e.gridToL1MedianM/1e3, gridRadial/1e3, gridClock/1e3, ...
    gridTransverse/1e3, e.gridCost, e.costMin);
end

function s = solveL1(obs, position, clockM)
% 与 epochL1Check 相同：从 LS 解出发的 IRLS，得到单历元 L1 最优。
x = [position(:); clockM];
for iter = 1:300
    d = x(1:3)' - obs.satPosition; r0 = vecnorm(d,2,2);
    residual = obs.range - r0 - x(4);
    w = 1./max(abs(residual),1e-4);
    G = [d./r0, ones(numel(r0),1)];
    step = (sqrt(w).*G) \ (sqrt(w).*residual);
    x = x + step;
    if norm(step) < 1e-4, break; end
end
s.position = x(1:3); s.clockM = x(4);
end

%% 批处理
function b = batchLandscape(cfg)
cfg.batch.methods = {'dpeGrid'};
assert(~cfg.batch.cn0Weight,'LuGRE:Config','costLandscape 只实现等权代价。');
res = runBatchPositioning(cfg);
obs = res.observations; epochs = res.gpsSeconds;
rotations = icrfToItrfRotation(epochs, cfg.frame);
[H, y] = linearize(res.dpeGrid.state, res.moon0, res.tRef, epochs, rotations, obs, cfg);
counts = accumarray(obs.epochIndex(:), 1, [numel(epochs) 1]);
search = counts(obs.epochIndex(:)) >= 2;
layout = epochLayout(obs.epochIndex(search));
T = whitening(H(search,:), layout);
HT = H(search,:)*T; y = y(search);
half = 20;
axisU = linspace(-half, half, 401);
[U1,U2] = meshgrid(axisU, axisU);
Ufull = zeros(6, numel(U1)); Ufull(1,:) = U1(:)'; Ufull(2,:) = U2(:)';
J = reshape(gridCost(HT, y, layout, Ufull, 2048), size(U1));
% 白化前：同一平面在物理坐标中的形状（u1、u2 方向对应 tRef 位置的偏移，单位 m）。
lengthM = [norm(T(1:3,1)), norm(T(1:3,2))];
b = struct('uAxis',axisU,'cost',J,'costMin',min(J(:)),'lattice',(-2:2)*half/2.5, ...
    'positionPerUnitM',lengthM,'observations',nnz(search),'epochs',nnz(counts>=2));
fprintf('%s 批处理：%d 条观测、%d 个多星历元；u1/u2 每单位对应 tRef 位置 %.1f / %.1f m\n', ...
    cfg.data.window, b.observations, b.epochs, lengthM);
end

% 以下三个函数与 src/batch/solveBatchGrid.m 中的局部函数相同（等权版本）。
function [H, y] = linearize(state, moon0, tRef, epochs, rotations, obs, cfg)
steps = [1e3*ones(3,1); ones(3,1)];
positionIcrf = propagateArc(state, moon0, tRef, epochs, cfg);
m = numel(epochs);
Phi = zeros(3,6,m);
for j = 1:6
    dx = zeros(6,1); dx(j) = steps(j);
    p = propagateArc(state + dx, moon0, tRef, epochs, cfg);
    Phi(:,j,:) = permute((p - positionIcrf)/steps(j), [2 3 1]);
end
positionEcef = zeros(m,3);
for k = 1:m, positionEcef(k,:) = (rotations(:,:,k)*positionIcrf(k,:)')'; end
idx = obs.epochIndex(:);
los = positionEcef(idx,:) - obs.satPosition;
distance = vecnorm(los,2,2);
unit = los./distance;
H = zeros(numel(idx),6);
for k = unique(idx)'
    rows = idx == k;
    H(rows,:) = unit(rows,:)*rotations(:,:,k)*Phi(:,:,k);
end
y = obs.range(:) - distance;
end

function T = whitening(H, layout)
wsum = accumarray(layout.group, 1);
Hc = H;
for j = 1:size(H,2)
    mu = accumarray(layout.group, H(:,j))./wsum;
    Hc(:,j) = H(:,j) - mu(layout.group);
end
n = vecnorm(Hc,2,1); n(n==0) = 1;
Hs = Hc./n;
[V, D] = eig((Hs'*Hs + (Hs'*Hs)')/2);
T = (V./n')./sqrt(diag(D)');
end

function cost = gridCost(HT, y, layout, U, chunkSize)
G = size(U,2); cost = zeros(1,G);
[mEff, K] = size(layout.order);
W = double(layout.order <= numel(y));
for first = 1:chunkSize:G
    cols = first:min(first+chunkSize-1, G);
    R = [y - HT*U(:,cols); nan(1,numel(cols))];
    Rp = reshape(R(layout.order,:), mEff, K, numel(cols));
    med = median(Rp, 2, 'omitnan');
    cost(cols) = reshape(sum(W.*abs(Rp - med), [1 2], 'omitnan'), 1, []);
end
end

function layout = epochLayout(epochIndex)
[~, ~, group] = unique(epochIndex);
counts = accumarray(group, 1);
mEff = numel(counts); K = max(counts);
[~, sorter] = sort(group);
order = zeros(mEff, K);
position = 1;
for k = 1:mEff
    order(k,1:counts(k)) = sorter(position:position+counts(k)-1)';
    position = position + counts(k);
end
order(order == 0) = numel(epochIndex) + 1;
layout = struct('group',group,'order',order);
end

%% 作图
function plotLandscape(land, outputFolder)
% 三栏：逐历元全景（km）、逐历元格点放大（m）、批处理白化平面。色标均为 ΔJ（m，对数）。
e = land.epoch; b = land.batch; name = strrep(land.window,'_0','');
font = 'PingFang SC';
fig = figure('Visible','off','Color','w','Units','pixels','Position',[0 0 1800 620]);
theme(fig,'light');
tl = tiledlayout(fig,1,3,'TileSpacing','loose','Padding','loose');
title(tl, sprintf('%s 真实观测的 L1 代价曲面', name),'FontName',font,'FontSize',18,'FontWeight','bold');
dot = [0.35 0.38 0.45]; star = [0.85 0.15 0.15]; gridColor = [1 0.55 0];

% 1. 逐历元全景
ax = nexttile(tl);
showCost(ax, e.radialM/1e3, e.clockM/1e3, e.cost - e.costMin);
hold(ax,'on');
h = plot(ax, 0, 0, 'p','MarkerSize',16,'MarkerFaceColor',star,'MarkerEdgeColor','w');
h(2) = plot(ax, e.gridRadialM/1e3, e.gridClockM/1e3, 'o','MarkerSize',10,'LineWidth',2.5,'Color',gridColor);
h(3) = plot(ax, e.lsRadialM/1e3, e.lsClockM/1e3, 'd','MarkerSize',8,'MarkerFaceColor','w','MarkerEdgeColor','k');
zr = e.zoomRadialM([1 end])/1e3; zc = e.zoomClockM([1 end])/1e3;
rectangle(ax,'Position',[zr(1) zc(1) diff(zr) diff(zc)],'EdgeColor','w','LineWidth',1.5,'LineStyle','--');
legend(ax, h, {'L1 最优','网格 DPE 解','LS 解'},'Location','northeast','FontName',font,'FontSize',10);
xlabel(ax,'径向位置偏移 (km)','FontName',font); ylabel(ax,'钟差偏移 c·dt (km)','FontName',font);
title(ax, sprintf('逐历元（第 %d 历元，%d 星）：径向-钟差平面', e.epochIndex, e.satCount), ...
    'FontName',font,'FontSize',13);
subtitle(ax, sprintf('网格解沿谷离 L1 最优 %.0f km；代价 %.0f m vs %.1f m', ...
    e.gridToL1M/1e3, e.gridCost, e.costMin),'FontName',font,'FontSize',10);

% 2. 放大：格点与谷底
ax = nexttile(tl);
showCost(ax, (e.zoomRadialM - e.gridRadialM)/1e3, (e.zoomClockM - e.gridClockM)/1e3, e.zoomCost - e.costMin);
hold(ax,'on');
[LR,LC] = meshgrid((e.latticeRadialM - e.gridRadialM)/1e3, (e.latticeClockM - e.gridClockM)/1e3);
inside = abs(LR) < 2.9*e.spacingM/1e3 & abs(LC) < 2.9*e.spacingM/1e3;
plot(ax, LR(inside), LC(inside), 'o','MarkerSize',7,'MarkerFaceColor','w','MarkerEdgeColor',dot);
plot(ax, 0, 0, 'o','MarkerSize',12,'LineWidth',2.5,'Color',gridColor);
xlabel(ax,'径向偏移 (km，相对网格解)','FontName',font); ylabel(ax,'钟差偏移 (km，相对网格解)','FontName',font);
% 格点处的代价增量：离谷轴几百米的代价差远大于沿谷几十 km 的代价差，最低格点因此沿谷乱跑。
latticeCost = interp2(e.zoomRadialM, e.zoomClockM, e.zoomCost - e.costMin, ...
    LR(inside)*1e3 + e.gridRadialM, LC(inside)*1e3 + e.gridClockM);
title(ax, sprintf('放大网格解附近：网格间距 %.1f km', e.spacingM/1e3),'FontName',font,'FontSize',13);
subtitle(ax, sprintf('格点处 ΔJ 中位数 %.0f m，远大于沿谷 %.0f km 的代价差 %.0f m', ...
    median(latticeCost), e.gridToL1M/1e3, e.gridCost - e.costMin),'FontName',font,'FontSize',10);

% 3. 批处理白化平面
ax = nexttile(tl);
showCost(ax, b.uAxis, b.uAxis, b.cost - b.costMin);
hold(ax,'on');
[L1,L2] = meshgrid(b.lattice);
plot(ax, L1(:), L2(:), 'o','MarkerSize',7,'MarkerFaceColor','w','MarkerEdgeColor',dot);
plot(ax, 0, 0, 'p','MarkerSize',16,'MarkerFaceColor',star,'MarkerEdgeColor','w');
xlabel(ax, sprintf('白化坐标 u_1（最弱方向，1 单位 ≈ %.0f m）', b.positionPerUnitM(1)),'FontName',font);
ylabel(ax, sprintf('白化坐标 u_2（1 单位 ≈ %.0f m）', b.positionPerUnitM(2)),'FontName',font);
title(ax,'批处理 dpeGrid：白化坐标','FontName',font,'FontSize',13);
subtitle(ax, sprintf('%d 条观测、%d 个历元；圆点为一级 5×5 网格；白化前此平面为 %.1f:1 的椭圆', ...
    b.observations, b.epochs, max(b.positionPerUnitM)/min(b.positionPerUnitM)), ...
    'FontName',font,'FontSize',10);

file = fullfile(outputFolder, [land.window '_landscape']);
exportgraphics(fig, [file '.png'], 'Resolution', 200);
exportgraphics(fig, [file '.pdf'], 'ContentType', 'vector');
close(fig);
end

function showCost(ax, x, y, dJ)
imagesc(ax, x, y, max(dJ, 0.1));
set(ax,'YDir','normal','ColorScale','log','FontName','PingFang SC','FontSize',11);
axis(ax,'equal'); axis(ax,'tight');
colormap(ax, flipud(sky(256)).^0.8);
cb = colorbar(ax); cb.Label.String = '\DeltaJ (m)'; cb.Label.FontSize = 12;
hold(ax,'on');
contour(ax, x, y, max(dJ,0.1), 10.^(0:0.5:5), 'EdgeColor',[0.2 0.2 0.3], 'EdgeAlpha',0.35, 'LineWidth',0.5);
end
