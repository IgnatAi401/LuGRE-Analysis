function check = epochL1Check(window, outputFolder)
% EPOCHL1CHECK 逐历元网格 DPE 的误差来源诊断：网格搜索是否找到 L1 代价的全局最优？
% 每个 NAV 历元比较三种解（同一组观测，轨道先验关闭）：
%   网格 DPE（solveEpochDpe）、精确 L1 最优（从 LS 解出发的 IRLS）、LS（L2 最优）。
% 若网格解代价≈精确 L1 代价，则误差来自 L1 代价本身（单历元只拟合 4 条观测），而非搜索失败。
% 例：epochL1Check('OP1_0')
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, genpath(fullfile(root,'src')));
if nargin < 2, outputFolder = fullfile(root,'results','method_compare'); end
cfg = lugreConfig(); cfg.data.window = window; cfg.enable.orbit = false;
if strcmp(window,'OP14_0'), cfg.run.skipEpochs = 3; end
[rawPath,navPath] = findTelemetryFiles(cfg.data.root,window);
nav = readTelemetryNav(navPath); nav = nav(cfg.run.skipEpochs+1:end);
raw = readTelemetryRaw(rawPath);
time = [nav.rxTime]'; n = numel(time);
[gps,galileo] = loadRinexFiles(selectRinexFiles(cfg,time));
raw = addTransmitTime(raw(ismember([raw.rxTime],time)),cfg.constants.c);
raw = addSatelliteStates(raw,gps,galileo,cfg.ephemeris.maxAgeS);
truth = loadTruthEcef(time, cfg.truth, cfg.frame);

blank = nan(n,1);
check = struct('window',window,'gpsSeconds',time,'satCount',blank,'pdop',blank, ...
    'gridCost',blank,'l1Cost',blank,'lsCost',blank,'gridError',blank,'l1Error',blank,'lsError',blank, ...
    'gridToL1M',blank,'l1ExactFits',blank,'gridRadial',blank,'l1Radial',blank,'lsRadial',blank);
dpePosition = [nav(1).posX; nav(1).posY; nav(1).posZ]; dpeClock = cfg.initial.clockBiasM;
lsPosition = dpePosition; lsClock = dpeClock;
for i = 1:n
    obs = buildObservations(selectEpochMeasurements(time(i),raw),dpePosition,cfg);
    check.satCount(i) = numel(obs.range);
    grid = solveEpochDpe(obs,dpePosition,dpeClock,nan(3,1),i,cfg);
    ls = solveEpochLs(obs,lsPosition,lsClock,cfg);
    if grid.valid, dpePosition = grid.position; dpeClock = grid.clockM; end
    if ls.valid, lsPosition = ls.position; lsClock = ls.clockM; end
    if ~(grid.valid && ls.valid), continue; end
    l1 = solveL1(obs, ls.position, ls.clockM);
    u = truth(i,:)'/norm(truth(i,:));
    cost = @(p,c) sum(abs(obs.range - vecnorm(p'-obs.satPosition,2,2) - c));
    check.gridCost(i) = cost(grid.position,grid.clockM);
    check.l1Cost(i) = cost(l1.position,l1.clockM);
    check.lsCost(i) = cost(ls.position,ls.clockM);
    check.gridError(i) = norm(grid.position-truth(i,:)'); check.gridRadial(i) = u'*(grid.position-truth(i,:)');
    check.l1Error(i) = norm(l1.position-truth(i,:)'); check.l1Radial(i) = u'*(l1.position-truth(i,:)');
    check.lsError(i) = norm(ls.position-truth(i,:)'); check.lsRadial(i) = u'*(ls.position-truth(i,:)');
    check.gridToL1M(i) = norm(grid.position-l1.position);
    check.l1ExactFits(i) = nnz(abs(l1.residual) < 0.5);
    G = [(ls.position'-obs.satPosition)./vecnorm(ls.position'-obs.satPosition,2,2), ones(numel(obs.range),1)];
    Q = inv(G'*G); check.pdop(i) = sqrt(trace(Q(1:3,1:3)));
end
good = isfinite(check.gridCost);
fprintf('%s：%d 历元。代价中位数 网格 %.1f / 精确L1 %.1f / LS %.1f m；网格解距 L1 最优中位数 %.0f m\n', ...
    window, nnz(good), median(check.gridCost(good)), median(check.l1Cost(good)), median(check.lsCost(good)), ...
    median(check.gridToL1M(good)));
fprintf('  真值误差中位数 网格 %.0f / 精确L1 %.0f / LS %.0f m；径向 RMS %.0f / %.0f / %.0f m；PDOP 中位数 %.0f；星数中位数 %d；L1 精确拟合观测数中位数 %d\n', ...
    median(check.gridError(good)), median(check.l1Error(good)), median(check.lsError(good)), ...
    rms(check.gridRadial(good)), rms(check.l1Radial(good)), rms(check.lsRadial(good)), ...
    median(check.pdop(good)), median(check.satCount(good)), median(check.l1ExactFits(good)));
if ~exist(outputFolder,'dir'), mkdir(outputFolder); end
save(fullfile(outputFolder,[window '_epochL1.mat']),'check');
end

function s = solveL1(obs, position, clockM)
% 单历元 L1 最优：IRLS（高斯-牛顿 + 1/|r| 权），从 LS 解出发；代价凸性弱，但 LS 解已在附近。
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
s.residual = obs.range - vecnorm(x(1:3)'-obs.satPosition,2,2) - x(4);
end
