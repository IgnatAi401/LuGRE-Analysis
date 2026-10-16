function summary = compareMethods(windows, outputFolder, options)
% COMPAREMETHODS 同一 OP 上运行全部定位方法并与真值对比，供汇报使用。
% 逐历元：网格 DPE（有/无轨道先验）、LS；批处理：dpeGrid、dpeIrls、ls；卡尔曼滤波：ekfOrbit、ekfCvBlend
% （LS 协方差初始化，另跑原项目固定协方差初始化作对照）；另有快照 LS、星上 NAV、轨道预测。
% 除各方法自身历元的统计外，另给出所有方法共同历元（NAV 历元 ∩ 真值有效）上的统计。
% 方法不可用时不运行，原因记入 info.status：没有 >=4 星历元的 OP 跳过全部方法；没有 NAV 文件时跳过逐历元方法；
% 没有 >=minInitSatellites 星历元时跳过滤波。逐历元方法按 NAV 连续 1 Hz 段分别处理（cfg.run.splitGaps）。
% 月面窗口（OP 编号 >= 38，着陆器静止在月面）：批处理自动用月固系静止模型（cfg.batch.dynamics='lunarStatic'），
% 滤波另加 ekfLunarStatic（月固系固定位置 + 钟差/钟漂）；其他 OP 的 ekfLunarStatic 记为 not applicable。
% options.skipBatchSolvers=true 时只算快照 LS，不运行 dpeGrid/dpeIrls/ls，状态记为 skipped
% （2026-10-13 曾用于月面长窗口：当时只有自由飞行模型，OP77/OP78 上均不收敛）。
% 例：compareMethods({'OP14_0','OP5_0'})
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, genpath(fullfile(root,'src')));
if nargin < 1, windows = {'OP14_0','OP5_0'}; end
if ischar(windows), windows = {windows}; end
if nargin < 2, outputFolder = fullfile(root,'results','method_compare'); end
if ~exist(outputFolder,'dir'), mkdir(outputFolder); end
if nargin < 3, options = struct(); end
if ~isfield(options,'skipBatchSolvers'), options.skipBatchSolvers = false; end
columns = {'window','method','count','medianM','rmseM','p95M', ...
    'commonCount','commonMedianM','commonRmseM','commonP95M','commonRadialRmseM','commonTransverseRmseM','seconds'};
methodNames = {'epochDpeOrbit','epochDpeFree','epochLs','orbitPrediction','batchDpeGrid','batchDpeIrls', ...
    'batchLs','ekfOrbit','ekfCvBlend','ekfOrbitFixedInit','ekfCvBlendFixedInit','ekfLunarStatic','snapshotLs','nav'};
rows = cell(0,numel(columns));
for w = 1:numel(windows)
    window = windows{w};
    cfg = lugreConfig(); cfg.data.window = window; cfg.run.plot = false; cfg.run.progressEvery = 0;
    cfg.run.splitGaps = true; % 例如 OP14 前 3 个 NAV 历元后有 178 s 空档，各段独立处理
    lunarSurface = sscanf(window,'OP%d') >= 38; % OP38 起着陆器静止在月面
    if lunarSurface, cfg.batch.dynamics = 'lunarStatic'; end
    windowStart = tic;
    [time, satCount, navPosition, hasNav] = surveyWindow(cfg);
    m = numel(time);
    status = struct();
    for k = 1:numel(methodNames), status.(methodNames{k}) = 'not run'; end
    status.nav = ternary(hasNav,'ok','unavailable: no NAV file');
    if ~lunarSurface, status.ekfLunarStatic = 'not applicable: lunar-surface model (OP>=38 only)'; end
    est = struct(); seconds = struct();
    for k = 1:numel(methodNames), seconds.(methodNames{k}) = NaN; end
    if hasNav, est.nav = navPosition; end
    count4 = nnz(satCount>=4);
    fprintf('\n===== %s：RAW %d 历元，>=4 星 %d，>=%d 星 %d，NAV %s =====\n', window, m, count4, ...
        cfg.filter.minInitSatellites, nnz(satCount>=cfg.filter.minInitSatellites), ternary(hasNav,'有','无'));
    epochOrbit = []; epochFree = []; batch = []; ekf = []; ekfFixed = [];
    if count4 == 0
        reason = 'unavailable: no epoch with >=4 satellites';
        for k = 1:numel(methodNames)-1, status.(methodNames{k}) = reason; end
    else
        % 逐历元（NAV 历元为轴，需要 NAV 初值）
        epochMethods = {'epochDpeOrbit','epochLs','orbitPrediction','epochDpeFree'};
        if ~hasNav
            for k = 1:numel(epochMethods), status.(epochMethods{k}) = 'unavailable: no NAV file'; end
        else
            [epochOrbit, secondsOrbit, message] = attempt(@()runEpochPositioning(cfg), window, '逐历元（轨道先验开）');
            if isempty(epochOrbit)
                [status.epochDpeOrbit, status.epochLs, status.orbitPrediction] = deal(message);
            else
                [status.epochDpeOrbit, status.epochLs, status.orbitPrediction] = deal('ok');
                seconds.epochDpeOrbit = secondsOrbit;
            end
            c = cfg; c.enable.orbit = false; c.enable.ls = false;
            [epochFree, secondsFree, message] = attempt(@()runEpochPositioning(c), window, '逐历元（轨道先验关）');
            status.epochDpeFree = message;
            if ~isempty(epochFree), status.epochDpeFree = 'ok'; seconds.epochDpeFree = secondsFree; end
        end
        % 批处理（初值来自快照 LS，需要 >=initMinEpochs 个 >=4 星历元）
        c = cfg;
        if options.skipBatchSolvers, c.batch.methods = {}; end
        [batch, secondsBatch, message] = attempt(@()runBatchPositioning(c), window, '批处理');
        if isempty(batch)
            [status.batchDpeGrid, status.batchDpeIrls, status.batchLs, status.snapshotLs] = deal(message);
        else
            assert(isequal(batch.gpsSeconds, time), '批处理与统计的 RAW 历元轴不一致。');
            status.snapshotLs = 'ok';
            est.snapshotLs = batch.init.snapshotPosition;
            batchNames = struct('dpeGrid','batchDpeGrid','dpeIrls','batchDpeIrls','ls','batchLs');
            for name = fieldnames(batchNames)'
                if ~isfield(batch, name{1})
                    status.(batchNames.(name{1})) = 'skipped: options.skipBatchSolvers';
                    continue;
                end
                s = batch.(name{1});
                status.(batchNames.(name{1})) = ternary(s.converged,'ok', ...
                    sprintf('ok (not converged after %d outer iterations)', s.iterations));
                est.(batchNames.(name{1})) = s.positionEcef;
            end
            if isfield(batch,'dpeGrid'), seconds.batchDpeGrid = batch.dpeGrid.search.seconds; end
        end
        % 卡尔曼滤波（需要 >=minInitSatellites 星的初始化历元）
        if nnz(satCount >= cfg.filter.minInitSatellites) == 0
            reason = sprintf('unavailable: no epoch with >=%d satellites', cfg.filter.minInitSatellites);
            [status.ekfOrbit, status.ekfCvBlend, status.ekfOrbitFixedInit, status.ekfCvBlendFixedInit] = deal(reason);
            if lunarSurface, status.ekfLunarStatic = reason; end
        else
            c = cfg;
            if lunarSurface, c.filter.methods{end+1} = 'ekfLunarStatic'; end
            [ekf, ~, message] = attempt(@()runFilterPositioning(c), window, '卡尔曼滤波');
            [status.ekfOrbit, status.ekfCvBlend] = deal(message);
            if lunarSurface, status.ekfLunarStatic = message; end
            if ~isempty(ekf)
                assert(isequal(ekf.gpsSeconds, time), '滤波与统计的 RAW 历元轴不一致。');
                [status.ekfOrbit, status.ekfCvBlend] = deal('ok');
                est.ekfOrbit = ekf.ekfOrbit.positionEcef; seconds.ekfOrbit = ekf.ekfOrbit.seconds;
                est.ekfCvBlend = ekf.ekfCvBlend.positionEcef; seconds.ekfCvBlend = ekf.ekfCvBlend.seconds;
                if lunarSurface
                    status.ekfLunarStatic = 'ok';
                    est.ekfLunarStatic = ekf.ekfLunarStatic.positionEcef;
                    seconds.ekfLunarStatic = ekf.ekfLunarStatic.seconds;
                end
            end
            c = cfg; c.filter.initCovariance = 'fixed';
            [ekfFixed, ~, message] = attempt(@()runFilterPositioning(c), window, '卡尔曼滤波（固定初始协方差）');
            [status.ekfOrbitFixedInit, status.ekfCvBlendFixedInit] = deal(message);
            if ~isempty(ekfFixed)
                [status.ekfOrbitFixedInit, status.ekfCvBlendFixedInit] = deal('ok');
                est.ekfOrbitFixedInit = ekfFixed.ekfOrbit.positionEcef;
                est.ekfCvBlendFixedInit = ekfFixed.ekfCvBlend.positionEcef;
            end
        end
    end

    % 统一到 RAW 历元轴，逐历元结果按 GPS 秒填入。
    inEpoch = false(m,1);
    if ~isempty(epochOrbit) || ~isempty(epochFree)
        epochAxis = epochOrbit; if isempty(epochAxis), epochAxis = epochFree; end
        [inEpoch, loc] = ismembertol(time, epochAxis.gpsSeconds, 1e-6, 'DataScale', 1);
        place = @(p) fillRows(nan(m,3), inEpoch, p(loc(inEpoch),:));
        if ~isempty(epochOrbit)
            est.epochDpeOrbit = place(maskInvalid(epochOrbit.dpe));
            est.epochLs = place(maskInvalid(epochOrbit.ls));
            est.orbitPrediction = place(epochOrbit.orbit.position);
        end
        if ~isempty(epochFree), est.epochDpeFree = place(maskInvalid(epochFree.dpe)); end
    elseif hasNav
        inEpoch = all(isfinite(navPosition),2);
    end
    truth = nan(m,3);
    if m > 0
        try
            [truthPosition, truthValid] = loadTruthEcef(time, cfg.truth, cfg.frame);
            truth(truthValid,:) = truthPosition(truthValid,:);
        catch ME
            warning('LuGRE:TruthUnavailable','%s 真值不可用：%s', window, ME.message);
        end
    end

    common = inEpoch & all(isfinite(truth),2);
    radialUnit = truth./vecnorm(truth,2,2);
    series = struct('gpsSeconds',time,'truth',truth,'satCount',satCount, ...
        'navEpoch',inEpoch,'common',common);
    for k = 1:numel(methodNames)
        name = methodNames{k};
        if ~isfield(est,name), continue; end
        e = est.(name) - truth;
        radial = sum(e.*radialUnit,2);
        transverse = vecnorm(e - radial.*radialUnit,2,2);
        d = vecnorm(e,2,2);
        series.(name) = struct('positionEcef',est.(name),'distanceM',d,'radialM',radial,'transverseM',transverse);
        own = isfinite(d); both = own & common;
        rows(end+1,:) = {window, name, nnz(own), safeStat(@median,d(own)), safeStat(@rms,d(own)), ...
            safeStat(@(x)prctile(x,95),d(own)), nnz(both), safeStat(@median,d(both)), safeStat(@rms,d(both)), ...
            safeStat(@(x)prctile(x,95),d(both)), safeStat(@rms,radial(both)), safeStat(@rms,transverse(both)), ...
            seconds.(name)}; %#ok<AGROW>
    end
    info = struct('window',window,'rawEpochs',m,'navEpochs',nnz(inEpoch),'hasNav',hasNav, ...
        'epochsWith4',count4,'satCountHistogram',accumarray(satCount+1,1)','status',status, ...
        'meanRangeRe',mean(vecnorm(truth,2,2),'omitnan')/cfg.constants.earthRadiusM, ...
        'batchDynamics',cfg.batch.dynamics,'secondsWindowTotal',NaN);
    if ~isempty(epochOrbit)
        info.epochSegments = nnz(epochOrbit.segmentStart);
        info.epochDpeBoundary = mean(epochOrbit.dpe.onBoundary(epochOrbit.dpe.valid));
        info.epochLsClockM = epochOrbit.ls.clockM;
    end
    if ~isempty(epochFree), info.epochDpeFreeBoundary = mean(epochFree.dpe.onBoundary(epochFree.dpe.valid)); end
    if ~isempty(batch)
        info.observations = numel(batch.observations.range);
        info.secondsBatchTotal = secondsBatch;
    end
    if ~isempty(batch) && all(isfield(batch,{'dpeGrid','dpeIrls','ls'}))
        info.gridPoints = batch.dpeGrid.search.totalPoints;
        info.gridLevels = sum([batch.dpeGrid.search.history.levels]);
        info.gridOuter = batch.dpeGrid.iterations;
        info.converged = struct('dpeGrid',batch.dpeGrid.converged,'dpeIrls',batch.dpeIrls.converged,'ls',batch.ls.converged);
        info.robustSigma = struct('dpeGrid',batch.dpeGrid.robustSigmaM,'dpeIrls',batch.dpeIrls.robustSigmaM, ...
            'ls',batch.ls.robustSigmaM);
        info.formalStdMedian = struct('dpeGrid',median(batch.dpeGrid.formalStd3dM), ...
            'dpeIrls',median(batch.dpeIrls.formalStd3dM),'ls',median(batch.ls.formalStd3dM));
        info.clockM = struct('dpeGrid',batch.dpeGrid.clockM,'dpeIrls',batch.dpeIrls.clockM,'ls',batch.ls.clockM);
        info.gridHistory = batch.dpeGrid.search.history;
        info.batchIterations = struct('dpeGrid',batch.dpeGrid.iterations,'dpeIrls',batch.dpeIrls.iterations, ...
            'ls',batch.ls.iterations);
        if lunarSurface % 月固系固定位置（m）：初值与三种批处理的解
            info.moonFixedState = struct('init',batch.init.state,'dpeGrid',batch.dpeGrid.state, ...
                'dpeIrls',batch.dpeIrls.state,'ls',batch.ls.state);
        end
    end
    if ~isempty(ekf)
        info.ekf = struct('initialEpoch',ekf.ekfOrbit.initialEpoch, ...
            'formalStd3dM',[ekf.ekfOrbit.formalStd3dM ekf.ekfCvBlend.formalStd3dM], ...
            'acceptedCount',[ekf.ekfOrbit.acceptedCount ekf.ekfCvBlend.acceptedCount], ...
            'rangeInnovationRmsM',[ekf.ekfOrbit.rangeInnovationRmsM ekf.ekfCvBlend.rangeInnovationRmsM], ...
            'cvBlendKalman',ekf.ekfCvBlend.kalmanPositionEcef,'cvBlendOrbit',ekf.ekfCvBlend.orbitPredictionEcef);
        if isfield(ekf,'ekfLunarStatic')
            s = ekf.ekfLunarStatic;
            info.ekfLunarStatic = struct('initialEpoch',s.initialEpoch,'formalStd3dM',s.formalStd3dM, ...
                'positionMoonFixed',s.positionMoonFixed,'acceptedCount',s.acceptedCount, ...
                'rangeInnovationRmsM',s.rangeInnovationRmsM);
        end
    end
    info.secondsWindowTotal = toc(windowStart);
    stats = cell2table(rows(strcmp(rows(:,1),window),:),'VariableNames',columns);
    statusTable = cell2table([fieldnames(status), struct2cell(status)],'VariableNames',{'method','status'});
    save(fullfile(outputFolder,[window '.mat']),'series','info','stats','statusTable');
    writetable(stats, fullfile(outputFolder,[window '_summary.csv']));
    writetable(statusTable, fullfile(outputFolder,[window '_status.csv']));
    fprintf('%s 完成，用时 %.0f s\n', window, info.secondsWindowTotal);
    disp(statusTable);
end
summary = cell2table(rows,'VariableNames',columns);
disp(summary);
end

function [time, satCount, navPosition, hasNav] = surveyWindow(cfg)
% 与批处理/滤波相同的 RAW 历元轴与每历元卫星数（每星一个信号、有星历），用于判断方法是否可用。
[rawPath, navPath] = findTelemetryFiles(cfg.data.root, cfg.data.window);
raw = readTelemetryRaw(rawPath);
raw = raw(isfinite([raw.prRaw]) & [raw.prRaw] > 0);
time = unique([raw.rxTime])';
raw = addTransmitTime(raw, cfg.constants.c);
[gps, galileo] = loadRinexFiles(selectRinexFiles(cfg, time));
raw = addSatelliteStates(raw, gps, galileo, cfg.ephemeris.maxAgeS);
obsCfg = cfg; obsCfg.enable.legacyIonosphere = false;
satCount = zeros(numel(time),1);
rawByEpoch = groupRawByEpoch(raw, time);
for k = 1:numel(time)
    satCount(k) = numel(buildObservations(selectEpochMeasurements(time(k), rawByEpoch{k}), zeros(3,1), obsCfg).range);
end
hasNav = ~isempty(navPath);
navPosition = nan(numel(time),3);
if hasNav
    nav = readTelemetryNav(navPath);
    [in, loc] = ismembertol(time, [nav.rxTime]', 1e-6, 'DataScale', 1);
    navPosition(in,:) = [[nav(loc(in)).posX]', [nav(loc(in)).posY]', [nav(loc(in)).posZ]'];
end
end

function [result, seconds, message] = attempt(run, window, label)
fprintf('\n===== %s %s =====\n', window, label);
t = tic;
try
    result = run(); message = 'ok';
catch ME
    result = []; message = ['failed: ' ME.identifier ' ' ME.message];
    warning('LuGRE:MethodFailed','%s %s 失败：%s', window, label, ME.message);
end
seconds = toc(t);
end

function value = safeStat(f, x)
if isempty(x), value = NaN; else, value = f(x); end
end

function p = maskInvalid(s)
p = s.position; p(~s.valid,:) = NaN;
end

function a = fillRows(a, mask, values)
a(mask,:) = values;
end

function out = ternary(condition, a, b)
if condition, out = a; else, out = b; end
end
