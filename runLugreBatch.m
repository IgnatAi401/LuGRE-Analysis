function results = runLugreBatch(cfg)
% RUNLUGREBATCH 多历元批处理 DPE 与同模型 LS 对照。
% 以 RAW 历元为时间轴，不足 4 颗星的历元也参与；初值来自快照 LS，不使用 NAV 位置/速度。
% NAV 若存在仅用于对比（及可选的初始化后备 cfg.batch.allowNavInit）。
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root,'function'), fullfile(root,'reader'));
validateBatchConfig(cfg);
[rawPath, navPath] = findWindowFiles(cfg.data.root, cfg.data.window);
raw = TxtParserRaw(rawPath);
raw = raw(isfinite([raw.prRaw]) & [raw.prRaw] > 0);
assert(~isempty(raw),'LuGRE:NoEpochs','RAW 没有有效伪距。');
epochs = unique([raw.rxTime])';
epochs = epochs(1:min(numel(epochs),cfg.run.maxEpochs));
raw = raw(ismember([raw.rxTime],epochs));
transmit = num2cell([raw.rxTime] - [raw.prRaw]/cfg.constants.c);
[raw.txTime] = transmit{:};
results.config = cfg;
results.gpsSeconds = epochs;
results.rinexFiles = selectRinexFiles(cfg, epochs);
[gps, galileo] = loadRinexFiles(results.rinexFiles);
raw = getPosRAW(raw, gps, galileo, cfg.ephemeris.maxAgeS);

% 每历元每颗卫星一个信号，模型与 main5 的 DPE/LS 相同（卫星钟差、Sagnac）。
m = numel(epochs);
obsCfg = cfg; obsCfg.enable.legacyIonosphere = false;
epochObs = cell(m,1);
for k = 1:m
    epochObs{k} = lugreObservations(getPos(epochs(k), raw), zeros(3,1), obsCfg);
end
counts = cellfun(@(o)numel(o.range), epochObs);
assert(any(counts > 0),'LuGRE:NoEpochs','没有可用观测（检查星历覆盖）。');
obs.epochIndex = repelem((1:m)', counts);
obs.satPosition = cell2mat(cellfun(@(o)o.position, epochObs, 'UniformOutput', false));
obs.range = cell2mat(cellfun(@(o)o.range, epochObs, 'UniformOutput', false));
obs.signalId = cell2mat(cellfun(@(o)o.signalId, epochObs, 'UniformOutput', false));
obs.svId = cell2mat(cellfun(@(o)o.svId, epochObs, 'UniformOutput', false));
obs.cn0 = cell2mat(cellfun(@(o)o.cn0, epochObs, 'UniformOutput', false));
obs.weight = ones(size(obs.range));
if cfg.batch.cn0Weight
    w = 10.^((obs.cn0-cfg.observation.cn0ReferenceDbHz)/10);
    w(~isfinite(w)) = cfg.observation.weightLimits(1);
    obs.weight = min(max(w,cfg.observation.weightLimits(1)),cfg.observation.weightLimits(2));
end
results.observationCount = counts;
results.observations = obs;
fprintf('%s: %d 个 RAW 历元，%d 个有观测，共 %d 条观测；>=4 星历元 %d 个。\n', ...
    cfg.data.window, m, nnz(counts>0), numel(obs.range), nnz(counts>=4));

observed = epochs(counts > 0);
[~,refIndex] = min(abs(epochs - (observed(1)+observed(end))/2));
tRef = epochs(refIndex);
rotations = zeros(3,3,m);
for k = 1:m, rotations(:,:,k) = lugreFrameRotation(epochs(k), cfg.frame); end

navPosition = [];
if ~isempty(navPath)
    nav = TxtParserNav(navPath);
    [in, loc] = ismembertol(epochs, [nav.rxTime]', 1e-6, 'DataScale', 1);
    navPosition = nan(m,3);
    navPosition(in,:) = [[nav(loc(in)).posX]', [nav(loc(in)).posY]', [nav(loc(in)).posZ]'];
end
results.navPosition = navPosition;
init = lugreBatchInit(epochObs, epochs, rotations, tRef, navPosition, cfg);
results.init = init;
fprintf('初值来源 %s：快照 LS 有效 %d 个历元。\n', init.source, nnz(all(isfinite(init.snapshotPosition),2)));
moon = getMoonEph(tRef, cfg);
moon0 = moon(2:7)'*1000; % JD 列不转换
results.tRef = tRef;
results.moon0 = moon0;

for k = 1:numel(cfg.batch.methods)
    method = cfg.batch.methods{k};
    solution = lugreBatchSolve(obs, epochs, rotations, init.state, moon0, tRef, method, cfg);
    results.(method) = solution;
    fprintf('%s 批处理: 收敛=%d, 外迭代 %d, 观测 %d, 参数 %d, 稳健残差尺度 %.2f m, 形式 3D std 中位数 %.0f m\n', ...
        upper(method), solution.converged, solution.iterations, solution.observationCount, ...
        solution.parameterCount, solution.robustSigmaM, median(solution.formalStd3dM));
    for j = 1:numel(solution.biasSignals)
        fprintf('    signal %d 相对 signal %d 偏差 %.2f m\n', solution.biasSignals(j), ...
            solution.referenceSignal, solution.signalBiasM(j));
    end
    if ~isempty(solution.grid)
        fprintf('    最弱方向网格最优偏移 %.0f m%s\n', solution.grid.bestOffsetM, ...
            ternary(solution.grid.onBoundary,'（在边界，需加大 gridHalfWidthM）',''));
    end
end

if cfg.truth.enabled
    try
        [truth, valid] = lugreTruthEcef(epochs, cfg.truth, cfg.frame);
        results.truth = evaluateAgainstTruth(results, truth, valid);
    catch ME
        warning('LuGRE:TruthUnavailable','真值评估跳过：%s', ME.message);
    end
end
if cfg.run.plot, plotLugreBatch(results); end
end

function truth = evaluateAgainstTruth(results, position, valid)
% 3D 误差并分解为地心径向（远距离几何的弱方向）与横向。
truth = struct('position',position,'valid',valid);
names = {}; estimates = {};
for method = results.config.batch.methods
    names{end+1} = ['batch_' method{1}]; estimates{end+1} = results.(method{1}).positionEcef; %#ok<AGROW>
end
names{end+1} = 'snapshotLs'; estimates{end+1} = results.init.snapshotPosition;
if ~isempty(results.navPosition)
    names{end+1} = 'nav'; estimates{end+1} = results.navPosition;
end
radialUnit = position./vecnorm(position,2,2);
fprintf('\n%-14s %6s %11s %11s %11s %13s %13s\n','真值对比','历元','中位数/m','RMSE/m','95%/m','径向RMSE/m','横向RMSE/m');
for k = 1:numel(names)
    e = estimates{k} - position;
    radial = sum(e.*radialUnit,2);
    transverse = vecnorm(e - radial.*radialUnit,2,2);
    d = vecnorm(e,2,2);
    good = isfinite(d);
    s = struct('errorM',e,'distanceM',d,'radialM',radial,'transverseM',transverse, ...
        'count',nnz(good),'medianM',median(d(good)),'rmseM',sqrt(mean(d(good).^2)), ...
        'p95M',prctile(d(good),95),'radialRmseM',sqrt(mean(radial(good).^2)), ...
        'transverseRmseM',sqrt(mean(transverse(good).^2)));
    truth.(names{k}) = s;
    fprintf('%-14s %6d %11.0f %11.0f %11.0f %13.0f %13.0f\n',names{k},s.count,s.medianM, ...
        s.rmseM,s.p95M,s.radialRmseM,s.transverseRmseM);
end
end

function [rawPath, navPath] = findWindowFiles(dataRoot, window)
folder = fullfile(dataRoot,'L0','TLM');
raw = dir(fullfile(folder,['TLM_RAW_*_' window '.txt']));
assert(isscalar(raw),'LuGRE:DataPath','窗口 %s 的 RAW 文件数为 %d（应为 1）。',window,numel(raw));
rawPath = fullfile(folder,raw.name);
nav = dir(fullfile(folder,['TLM_NAV_*_' window '.txt']));
navPath = '';
if isscalar(nav), navPath = fullfile(folder,nav.name); end
end

function validateBatchConfig(cfg)
b = cfg.batch;
assert(iscellstr(b.methods) && ~isempty(b.methods) && all(ismember(b.methods,{'dpe','ls'})), ...
    'LuGRE:Config','batch.methods 只能包含 dpe、ls。');
assert(any(strcmp(b.clockModel,{'polynomial','spline'})),'LuGRE:Config', ...
    'batch.clockModel 必须为 polynomial 或 spline。');
validateattributes(b.clockDegree,{'double'},{'scalar','integer','nonnegative'});
validateattributes(b.clockKnotSpacingS,{'double'},{'scalar','positive','finite'});
validateattributes(b.gridHalfWidthM,{'double'},{'scalar','positive','finite'});
validateattributes(b.gridStepM,{'double'},{'scalar','positive','finite'});
validateattributes(b.maxOuterIterations,{'double'},{'scalar','integer','positive'});
validateattributes(b.irlsEpsilonM,{'double'},{'scalar','positive','finite'});
validateattributes(cfg.run.maxEpochs,{'double'},{'scalar','positive'});
end

function out = ternary(condition, a, b)
if condition, out = a; else, out = b; end
end
