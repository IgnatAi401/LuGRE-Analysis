function summary = compareMethods(windows, outputFolder)
% COMPAREMETHODS 同一 OP 上运行全部定位方法并与真值对比，供汇报使用。
% 逐历元：网格 DPE（有/无轨道先验）、LS；批处理：dpeGrid、dpeIrls、ls；卡尔曼滤波：ekfOrbit、ekfCvBlend
% （LS 协方差初始化，另跑原项目固定协方差初始化作对照）；另有快照 LS、星上 NAV、轨道预测。
% 除各方法自身历元的统计外，另给出所有方法共同历元（NAV 历元 ∩ 真值有效）上的统计。
% 例：compareMethods({'OP14_0','OP5_0'})
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, genpath(fullfile(root,'src')));
if nargin < 1, windows = {'OP14_0','OP5_0'}; end
if nargin < 2, outputFolder = fullfile(root,'results','method_compare'); end
if ~exist(outputFolder,'dir'), mkdir(outputFolder); end
columns = {'window','method','count','medianM','rmseM','p95M', ...
    'commonCount','commonMedianM','commonRmseM','commonP95M','commonRadialRmseM','commonTransverseRmseM','seconds'};
rows = {};
for w = 1:numel(windows)
    window = windows{w};
    cfg = lugreConfig(); cfg.data.window = window; cfg.run.plot = false; cfg.run.progressEvery = 0;

    fprintf('\n===== %s 逐历元（轨道先验开）=====\n', window);
    if strcmp(window,'OP14_0'), cfg.run.skipEpochs = 3; end % 前 3 个 NAV 历元后有 178 s 空档
    t = tic; epochOrbit = runEpochPositioning(cfg); secondsOrbit = toc(t);
    c = cfg; c.enable.orbit = false; c.enable.ls = false;
    fprintf('\n===== %s 逐历元（轨道先验关）=====\n', window);
    t = tic; epochFree = runEpochPositioning(c); secondsFree = toc(t);
    fprintf('\n===== %s 批处理 =====\n', window);
    t = tic; batch = runBatchPositioning(cfg); secondsBatch = toc(t);
    fprintf('\n===== %s 卡尔曼滤波 =====\n', window);
    ekf = runFilterPositioning(cfg);
    c = cfg; c.filter.initCovariance = 'fixed';
    ekfFixed = runFilterPositioning(c);
    assert(isequal(ekf.gpsSeconds, batch.gpsSeconds), '滤波与批处理的 RAW 历元轴不一致。');

    % 统一到批处理的 RAW 历元轴，逐历元结果按 GPS 秒填入。
    time = batch.gpsSeconds; m = numel(time);
    truth = batch.truth.position;
    [inEpoch, loc] = ismembertol(time, epochOrbit.gpsSeconds, 1e-6, 'DataScale', 1);
    place = @(p) fillRows(nan(m,3), inEpoch, p(loc(inEpoch),:));
    est = struct();
    est.epochDpeOrbit = place(maskInvalid(epochOrbit.dpe));
    est.epochDpeFree = place(maskInvalid(epochFree.dpe));
    est.epochLs = place(maskInvalid(epochOrbit.ls));
    est.orbitPrediction = place(epochOrbit.orbit.position);
    est.batchDpeGrid = batch.dpeGrid.positionEcef;
    est.batchDpeIrls = batch.dpeIrls.positionEcef;
    est.batchLs = batch.ls.positionEcef;
    est.ekfOrbit = ekf.ekfOrbit.positionEcef;
    est.ekfCvBlend = ekf.ekfCvBlend.positionEcef;
    est.ekfOrbitFixedInit = ekfFixed.ekfOrbit.positionEcef;
    est.ekfCvBlendFixedInit = ekfFixed.ekfCvBlend.positionEcef;
    est.snapshotLs = batch.init.snapshotPosition;
    est.nav = batch.navPosition;
    seconds = struct('epochDpeOrbit',secondsOrbit,'epochDpeFree',secondsFree,'epochLs',NaN, ...
        'orbitPrediction',NaN,'batchDpeGrid',batch.dpeGrid.search.seconds,'batchDpeIrls',NaN, ...
        'batchLs',NaN,'ekfOrbit',ekf.ekfOrbit.seconds,'ekfCvBlend',ekf.ekfCvBlend.seconds, ...
        'ekfOrbitFixedInit',NaN,'ekfCvBlendFixedInit',NaN,'snapshotLs',NaN,'nav',NaN);

    common = inEpoch & all(isfinite(truth),2);
    radialUnit = truth./vecnorm(truth,2,2);
    names = fieldnames(est);
    series = struct('gpsSeconds',time,'truth',truth,'satCount',batch.observationCount, ...
        'navEpoch',inEpoch,'common',common);
    for k = 1:numel(names)
        e = est.(names{k}) - truth;
        radial = sum(e.*radialUnit,2);
        transverse = vecnorm(e - radial.*radialUnit,2,2);
        d = vecnorm(e,2,2);
        series.(names{k}) = struct('distanceM',d,'radialM',radial,'transverseM',transverse);
        own = isfinite(d); both = own & common;
        rows(end+1,:) = {window, names{k}, nnz(own), median(d(own)), rms(d(own)), prctile(d(own),95), ...
            nnz(both), median(d(both)), rms(d(both)), prctile(d(both),95), rms(radial(both)), ...
            rms(transverse(both)), seconds.(names{k})}; %#ok<AGROW>
    end
    info = struct('window',window,'rawEpochs',m,'navEpochs',nnz(inEpoch), ...
        'epochsWith4',nnz(batch.observationCount>=4),'observations',numel(batch.observations.range), ...
        'meanRangeRe',mean(vecnorm(truth,2,2),'omitnan')/cfg.constants.earthRadiusM, ...
        'epochDpeBoundary',mean(epochOrbit.dpe.onBoundary(epochOrbit.dpe.valid)), ...
        'epochDpeFreeBoundary',mean(epochFree.dpe.onBoundary(epochFree.dpe.valid)), ...
        'gridPoints',batch.dpeGrid.search.totalPoints,'gridLevels',sum([batch.dpeGrid.search.history.levels]), ...
        'gridOuter',batch.dpeGrid.iterations,'secondsBatchTotal',secondsBatch, ...
        'robustSigma',struct('dpeGrid',batch.dpeGrid.robustSigmaM,'dpeIrls',batch.dpeIrls.robustSigmaM, ...
        'ls',batch.ls.robustSigmaM), ...
        'formalStdMedian',struct('dpeGrid',median(batch.dpeGrid.formalStd3dM), ...
        'dpeIrls',median(batch.dpeIrls.formalStd3dM),'ls',median(batch.ls.formalStd3dM)), ...
        'clockM',struct('dpeGrid',batch.dpeGrid.clockM,'dpeIrls',batch.dpeIrls.clockM,'ls',batch.ls.clockM), ...
        'gridHistory',batch.dpeGrid.search.history, ...
        'ekf',struct('initialEpoch',ekf.ekfOrbit.initialEpoch,'formalStd3dM',[ekf.ekfOrbit.formalStd3dM ekf.ekfCvBlend.formalStd3dM], ...
        'acceptedCount',[ekf.ekfOrbit.acceptedCount ekf.ekfCvBlend.acceptedCount], ...
        'rangeInnovationRmsM',[ekf.ekfOrbit.rangeInnovationRmsM ekf.ekfCvBlend.rangeInnovationRmsM], ...
        'cvBlendKalman',ekf.ekfCvBlend.kalmanPositionEcef,'cvBlendOrbit',ekf.ekfCvBlend.orbitPredictionEcef));
    stats = cell2table(rows(strcmp(rows(:,1),window),:),'VariableNames',columns);
    save(fullfile(outputFolder,[window '.mat']),'series','info','stats');
    writetable(stats, fullfile(outputFolder,[window '_summary.csv']));
end
summary = cell2table(rows,'VariableNames',columns);
disp(summary);
end

function p = maskInvalid(s)
p = s.position; p(~s.valid,:) = NaN;
end

function a = fillRows(a, mask, values)
a(mask,:) = values;
end
