function results = runLugreAnalysis(cfg)
% RUNLUGREANALYSIS 独立执行 DPE/LS，保留 NAV 辅助的分段动力学方案。
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root,'function'), fullfile(root,'reader'));
validateConfig(cfg);
[~,rawPath,navPath] = FindTelemetryTxtFiles(cfg.data.root,cfg.data.window);
assert(~isempty(rawPath) && isfile(rawPath) && ~isempty(navPath) && isfile(navPath), ...
    'LuGRE:DataPath','未找到窗口 %s 的 RAW/NAV 文件。',cfg.data.window);
nav = TxtParserNav(navPath);
raw = TxtParserRaw(rawPath);
n = min(numel(nav),cfg.run.maxEpochs);
assert(n > 0, 'LuGRE:NoEpochs','NAV 没有历元。');
nav = nav(1:n);
time = [nav.rxTime]';
assert(all(isfinite(time)) && all(diff(time)>0),'LuGRE:EpochOrder','NAV 时间必须严格递增。');
assert(all(abs(diff(time)-1/cfg.data.sampleRateHz)<1e-6), ...
    'LuGRE:Sampling','NAV 时间间隔不符合配置的 1 Hz 数据约定。');
[gps,galileo] = read_gnss_rinex(cfg.data.rinex);
raw = getPosRAW(sortRAW(nav,raw),gps,galileo);
reference = [[nav.posX]',[nav.posY]',[nav.posZ]'];
results.config = cfg;
results.nav = nav;
results.gpsSeconds = time;
results.referencePosition = reference;
blank = struct('position',nan(n,3),'clockM',nan(n,1), 'valid',false(n,1));
results.dpe = blank; results.ls = blank;
results.dpe.onBoundary = false(n,1);
results.dpe.cost = nan(n,1);
results.ls.iterations = zeros(n,1);
results.observationCount = zeros(n,1);
results.orbit.position = nan(n,3);
results.orbit.enabled = cfg.enable.orbit;
results.moon = [];
if cfg.enable.orbit
    % 先验证工具箱/EOP，再发网络请求；不静默退回旧 GMST 近似。
    lugreFrameRotation(time(1),cfg.frame);
    moon = getMoonEph(time,cfg);
    moon(:,2:7) = moon(:,2:7)*1000; % JD_UTC 保持不变
    results.moon = moon;
end
origin = reference(1,:)' + cfg.initial.positionOffsetM(:);
dpePosition = origin; lsPosition = origin;
dpeClock = cfg.initial.clockBiasM; lsClock = dpeClock;
anchorState = []; anchorTime = NaN; lastReset = -Inf;
absTol = repmat([repmat(cfg.orbit.positionAbsToleranceM,3,1); ...
    repmat(cfg.orbit.velocityAbsToleranceMps,3,1)],2,1);
odeOptions = odeset('RelTol',cfg.orbit.relativeTolerance,'AbsTol',absTol);
for i = 1:n
    predicted = nan(3,1);
    if cfg.enable.orbit
        if isempty(anchorState)
            predicted = origin;
        else
            [~,states] = ode45(@(t,x)lugreDynamics(t,x,anchorTime,cfg), ...
                [0,time(i)-anchorTime],anchorState,odeOptions);
            [predicted,~] = j2000_to_ecef(states(end,1:3)',states(end,4:6)',time(i),cfg.frame);
        end
        results.orbit.position(i,:) = predicted';
    end
    % 两种方法使用完全相同的有效观测及传播修正。
    obs = lugreObservations(getPos(time(i),raw),origin,cfg);
    results.observationCount(i) = numel(obs.range);
    if cfg.enable.dpe
        prior = predicted;
        if i == 1, prior(:) = NaN; end % 第一历元与旧算法一样，不加先验
        solution = lugreDpe(obs,dpePosition,dpeClock,prior,i,cfg);
        results.dpe.position(i,:) = solution.position';
        results.dpe.clockM(i) = solution.clockM;
        results.dpe.valid(i) = solution.valid;
        results.dpe.onBoundary(i) = solution.onBoundary;
        results.dpe.cost(i) = solution.cost;
        if solution.valid
            dpePosition = solution.position; dpeClock = solution.clockM;
        end
    end
    if cfg.enable.ls
        solution = lugreLs(obs,lsPosition,lsClock,cfg);
        results.ls.position(i,:) = solution.position';
        results.ls.clockM(i) = solution.clockM;
        results.ls.valid(i) = solution.valid;
        results.ls.iterations(i) = solution.iterations;
        if solution.valid
            lsPosition = solution.position; lsClock = solution.clockM;
        end
    end
    if cfg.enable.dpe
        origin = dpePosition;
        feedbackValid = results.dpe.valid(i);
    else
        origin = lsPosition;
        feedbackValid = results.ls.valid(i);
    end
    % 双方法同时开启时动力学反馈取 DPE；LS-only 模式取 LS。
    % 保留用户指定的 NAV 速度辅助；不将预测误称作独立真值。
    if cfg.enable.orbit && (isempty(anchorState) || ...
            (time(i)-lastReset >= cfg.orbit.resetIntervalS && feedbackValid))
        velocity = [nav(i).velX;nav(i).velY;nav(i).velZ];
        [r,v] = ecef_to_j2000(origin,velocity,time(i),cfg.frame);
        anchorState = [r;v;moon(i,2:7)'];
        anchorTime = time(i); lastReset = time(i);
    end
    if cfg.run.progressEvery > 0 && (mod(i,cfg.run.progressEvery)==0 || i==n)
        fprintf('Processing %d / %d\n',i,n);
    end
end
for method = {'dpe','ls'}
    name = method{1};
    results.(name).differenceNavM = vecnorm(results.(name).position-reference,2,2);
    if cfg.enable.(name)
        values = results.(name).differenceNavM(results.(name).valid);
        fprintf('%s: 有效 %d/%d, 与 NAV 距离差中位数 %.3f m\n', ...
            upper(name),numel(values),n,median(values));
    end
end
results.orbit.differenceNavM = vecnorm(results.orbit.position-reference,2,2);
if isfield(cfg,'truth') && cfg.truth.enabled
    results = compareLugreTruth(results,cfg.truth);
end
if cfg.run.plot, plotLugreResults(results); end
end

function validateConfig(cfg)
assert(cfg.enable.dpe || cfg.enable.ls,'LuGRE:NoMethod','DPE 和 LS 至少开启一种。');
validateattributes(cfg.run.maxEpochs,{'double'},{'scalar','positive'});
assert(isinf(cfg.run.maxEpochs) || fix(cfg.run.maxEpochs)==cfg.run.maxEpochs, ...
    'LuGRE:Config','maxEpochs 必须为正整数或 Inf。');
assert(cfg.data.sampleRateHz==1,'LuGRE:Config','此数据集采样率固定为 1 Hz。');
if ~isempty(cfg.dpe.refineStepsM)
    validateattributes(cfg.dpe.refineStepsM,{'double'},{'row','positive','finite'});
end
assert(all(diff(cfg.dpe.refineStepsM)<0),'LuGRE:Config','精化步长必须递减。');
validateattributes(cfg.dpe.rankTolerance,{'double'},{'scalar','positive','finite'});
validateattributes(cfg.dpe.chunkSize,{'double'},{'scalar','integer','positive'});
validateattributes(cfg.orbit.resetIntervalS,{'double'},{'scalar','positive','finite'});
validateattributes(cfg.frame.derivativeStepS,{'double'},{'scalar','positive','finite'});
validateattributes(cfg.moon.stepMinutes,{'double'},{'scalar','integer','positive'});
validateattributes(cfg.dpe.orbitWeight,{'double'},{'scalar','nonnegative','finite'});
validateattributes(cfg.ls.maxIterations,{'double'},{'scalar','integer','positive'});
for key = {'coarseOffsetsM','fineOffsetsM','clockOffsetsM'}
    values = cfg.dpe.(key{1});
    validateattributes(values,{'double'},{'row','nonempty','finite'});
    assert(all(diff(values)>0),'LuGRE:Config','网格偏移必须严格递增。');
end
end
