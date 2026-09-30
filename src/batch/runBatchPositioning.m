function results = runBatchPositioning(cfg)
% RUNBATCHPOSITIONING 多历元批处理：dpeIrls、dpeGrid 与同模型 LS 对照。
% 以 RAW 历元为时间轴，不足 4 颗星的历元也参与；初值来自快照 LS，不使用 NAV 位置/速度。
% NAV 若存在仅用于对比（及可选的初始化后备 cfg.batch.allowNavInit）。
% dpeGrid 用 solveBatchGrid（每历元钟差），dpeIrls 与 ls 用 solveBatchIrls（样条/多项式钟差）。
validateBatchConfig(cfg);
[rawPath, navPath] = findTelemetryFiles(cfg.data.root, cfg.data.window);
raw = readTelemetryRaw(rawPath);
raw = raw(isfinite([raw.prRaw]) & [raw.prRaw] > 0);
assert(~isempty(raw),'LuGRE:NoEpochs','RAW 没有有效伪距。');
epochs = unique([raw.rxTime])';
epochs = epochs(1:min(numel(epochs),cfg.run.maxEpochs));
raw = raw(ismember([raw.rxTime],epochs));
raw = addTransmitTime(raw, cfg.constants.c);
results.config = cfg;
results.gpsSeconds = epochs;
results.rinexFiles = selectRinexFiles(cfg, epochs);
[gps, galileo] = loadRinexFiles(results.rinexFiles);
raw = addSatelliteStates(raw, gps, galileo, cfg.ephemeris.maxAgeS);

% 每历元每颗卫星一个信号，观测模型与逐历元 DPE/LS 相同（卫星钟差、Sagnac）。
m = numel(epochs);
obsCfg = cfg; obsCfg.enable.legacyIonosphere = false;
epochObs = cell(m,1);
for k = 1:m
    epochObs{k} = buildObservations(selectEpochMeasurements(epochs(k), raw), zeros(3,1), obsCfg);
end
counts = cellfun(@(o)numel(o.range), epochObs);
assert(any(counts > 0),'LuGRE:NoEpochs','没有可用观测（检查星历覆盖）。');
obs.epochIndex = repelem((1:m)', counts);
obs.satPosition = cell2mat(cellfun(@(o)o.satPosition, epochObs, 'UniformOutput', false));
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
rotations = icrfToItrfRotation(epochs, cfg.frame);

navPosition = [];
if ~isempty(navPath)
    nav = readTelemetryNav(navPath);
    [in, loc] = ismembertol(epochs, [nav.rxTime]', 1e-6, 'DataScale', 1);
    navPosition = nan(m,3);
    navPosition(in,:) = [[nav(loc(in)).posX]', [nav(loc(in)).posY]', [nav(loc(in)).posZ]'];
end
results.navPosition = navPosition;
init = initBatchState(epochObs, epochs, rotations, tRef, navPosition, cfg);
results.init = init;
fprintf('初值来源 %s：快照 LS 有效 %d 个历元。\n', init.source, nnz(all(isfinite(init.snapshotPosition),2)));
moon = fetchMoonEphemeris(tRef, cfg);
moon0 = moon(2:7)'*1000; % JD 列不转换
results.tRef = tRef;
results.moon0 = moon0;

for k = 1:numel(cfg.batch.methods)
    method = cfg.batch.methods{k};
    if strcmp(method,'dpeGrid')
        solution = solveBatchGrid(obs, epochs, rotations, init.state, moon0, tRef, cfg);
    else
        solution = solveBatchIrls(obs, epochs, rotations, init.state, moon0, tRef, method, cfg);
    end
    results.(method) = solution;
    fprintf('%s 批处理: 收敛=%d, 外迭代 %d, 观测 %d, 参数 %d, 稳健残差尺度 %.2f m, 形式 3D std 中位数 %.0f m\n', ...
        upper(method), solution.converged, solution.iterations, solution.observationCount, ...
        solution.parameterCount, solution.robustSigmaM, median(solution.formalStd3dM));
    for j = 1:numel(solution.biasSignals)
        fprintf('    signal %d 相对 signal %d 偏差 %.2f m\n', solution.biasSignals(j), ...
            solution.referenceSignal, solution.signalBiasM(j));
    end
    if isfield(solution,'search')
        h = solution.search.history;
        fprintf('    网格：%d 次外迭代共 %d 层、%.2e 点，边界移动 %d 次，耗时 %.1f s\n', numel(h), ...
            sum([h.levels]), solution.search.totalPoints, sum([h.boundaryMoves]), solution.search.seconds);
    end
    if ~isempty(solution.weakSearch)
        fprintf('    最弱方向搜索最优偏移 %.0f m%s\n', solution.weakSearch.bestOffsetM, ...
            ternary(solution.weakSearch.onBoundary,'（在边界，需加大 weakSearch.halfWidthM）',''));
    end
end

if cfg.truth.enabled
    try
        results.truth = compareBatchTruth(results);
    catch ME
        warning('LuGRE:TruthUnavailable','真值评估跳过：%s', ME.message);
    end
end
if cfg.run.plot, plotBatchResults(results); end
end

function validateBatchConfig(cfg)
b = cfg.batch;
assert(iscellstr(b.methods) && ~isempty(b.methods) && all(ismember(b.methods,{'dpeIrls','ls','dpeGrid'})), ...
    'LuGRE:Config','batch.methods 只能包含 dpeIrls、ls、dpeGrid。');
assert(any(strcmp(b.clockModel,{'polynomial','spline'})),'LuGRE:Config', ...
    'batch.clockModel 必须为 polynomial 或 spline。');
validateattributes(b.clockDegree,{'double'},{'scalar','integer','nonnegative'});
validateattributes(b.clockKnotSpacingS,{'double'},{'scalar','positive','finite'});
validateattributes(b.weakSearch.outerIterations,{'double'},{'scalar','integer','nonnegative'});
validateattributes(b.weakSearch.halfWidthM,{'double'},{'scalar','positive','finite'});
validateattributes(b.weakSearch.stepM,{'double'},{'scalar','positive','finite'});
validateattributes(b.maxOuterIterations,{'double'},{'scalar','integer','positive'});
validateattributes(b.irlsEpsilonM,{'double'},{'scalar','positive','finite'});
validateattributes(cfg.run.maxEpochs,{'double'},{'scalar','positive'});
end

function out = ternary(condition, a, b)
if condition, out = a; else, out = b; end
end
