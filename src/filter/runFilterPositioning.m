function results = runFilterPositioning(cfg)
% RUNFILTERPOSITIONING 扩展卡尔曼滤波定位（ekfOrbit、ekfCvBlend，见 solveFilterEkf）。
% 时间轴与批处理相同：全部 RAW 历元，1~3 星的历元也做量测更新；不使用 NAV 位置/速度，NAV 仅用于对比。
validateFilterConfig(cfg);
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

m = numel(epochs);
meas = cell(m,1);
for k = 1:m
    meas{k} = buildFilterMeasurements(selectEpochMeasurements(epochs(k), raw), cfg);
end
counts = cellfun(@(o)numel(o.range), meas);
results.observationCount = counts;
fprintf('%s: %d 个 RAW 历元，%d 个有观测；>=4 星历元 %d 个。\n', cfg.data.window, m, nnz(counts>0), nnz(counts>=4));

[rotations, rates] = icrfToItrfRotation(epochs, cfg.frame);
step = cfg.moon.stepMinutes*60;
grid = (floor(epochs(1)/step)*step - step : step : ceil(epochs(end)/step)*step + step)';
moonTable = fetchMoonEphemeris(grid, cfg);
moon = struct('gpsSeconds',grid,'state',moonTable(:,2:7)*1000); % JD 列不转换

results.navPosition = [];
if ~isempty(navPath)
    nav = readTelemetryNav(navPath);
    [in, loc] = ismembertol(epochs, [nav.rxTime]', 1e-6, 'DataScale', 1);
    results.navPosition = nan(m,3);
    results.navPosition(in,:) = [[nav(loc(in)).posX]', [nav(loc(in)).posY]', [nav(loc(in)).posZ]'];
end

for k = 1:numel(cfg.filter.methods)
    method = cfg.filter.methods{k};
    t = tic;
    solution = solveFilterEkf(meas, epochs, rotations, rates, moon, method, cfg);
    solution.seconds = toc(t);
    results.(method) = solution;
    fprintf('%s: 第 %d 个历元初始化，量测接受率 %.1f%%，伪距新息 RMS 中位数 %.1f m，形式 3D std 中位数 %.0f m，%.1f s\n', ...
        upper(method), solution.initialEpoch, 100*sum(solution.acceptedCount)/max(1,2*sum(counts(solution.initialEpoch:end))), ...
        median(solution.rangeInnovationRmsM,'omitnan'), median(solution.formalStd3dM,'omitnan'), solution.seconds);
end

if cfg.truth.enabled
    try
        [position, valid] = loadTruthEcef(epochs, cfg.truth, cfg.frame);
        names = cfg.filter.methods;
        estimates = cellfun(@(n)results.(n).positionEcef, names, 'UniformOutput', false);
        if ~isempty(results.navPosition)
            names{end+1} = 'nav'; estimates{end+1} = results.navPosition;
        end
        results.truth = summarizeTruthErrors(names, estimates, position);
        results.truth.valid = valid;
    catch ME
        warning('LuGRE:TruthUnavailable','真值评估跳过：%s', ME.message);
    end
end
if cfg.run.plot, plotFilterResults(results); end
end

function validateFilterConfig(cfg)
f = cfg.filter;
assert(iscellstr(f.methods) && ~isempty(f.methods) && all(ismember(f.methods,{'ekfOrbit','ekfCvBlend'})), ...
    'LuGRE:Config','filter.methods 只能包含 ekfOrbit、ekfCvBlend。');
validateattributes(f.minInitSatellites,{'double'},{'scalar','integer','>=',4});
validateattributes(f.initStd,{'double'},{'numel',8,'positive','finite'});
assert(any(strcmp(f.initCovariance,{'fixed','ls'})),'LuGRE:Config','filter.initCovariance 必须为 fixed 或 ls。');
validateattributes(f.ekfCvBlend.kalmanWeight,{'double'},{'scalar','>=',0,'<=',1});
validateattributes(f.maxStepS,{'double'},{'scalar','positive','finite'});
validateattributes(cfg.run.maxEpochs,{'double'},{'scalar','positive'});
end
