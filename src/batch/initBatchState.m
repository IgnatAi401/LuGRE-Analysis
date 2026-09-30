function init = initBatchState(epochObs, epochs, rotations, tRef, navPosition, cfg)
% INITBATCHSTATE 批处理初值：各历元快照 LS（自地心起算、不用 NAV），
% 转到 ICRF 后对时间做稳健（L1）线性拟合，得到 tRef 的 [r0;v0]。
% navPosition：与 epochs 对齐的 NAV 位置（无则为空），仅在 allowNavInit 时作后备。
opts = cfg.batch;
m = numel(epochs);
snapshot = nan(m,3); snapshotClock = nan(m,1);
lsCfg = cfg; lsCfg.enable.cn0Ls = false;
lsCfg.observation.minCount = max(4,opts.initMinSatellites);
previous = zeros(3,1); previousClock = 0;
for k = 1:m
    obs = epochObs{k};
    if numel(obs.range) < lsCfg.observation.minCount, continue; end
    s = solveEpochLs(obs, previous, previousClock, lsCfg);
    if ~s.valid && any(previous), s = solveEpochLs(obs, zeros(3,1), 0, lsCfg); end
    if s.valid
        snapshot(k,:) = s.position'; snapshotClock(k) = s.clockM;
        previous = s.position; previousClock = s.clockM;
    end
end
source = 'snapshotLs';
fitPosition = snapshot;
if nnz(all(isfinite(fitPosition),2)) < opts.initMinEpochs
    assert(opts.allowNavInit && ~isempty(navPosition) && ...
        nnz(all(isfinite(navPosition),2)) >= opts.initMinEpochs, 'LuGRE:BatchInit', ...
        '仅 %d 个历元有快照解（至少 %d 颗星），不足以初始化；可开启 cfg.batch.allowNavInit。', ...
        nnz(all(isfinite(snapshot),2)), opts.initMinSatellites);
    fitPosition = navPosition; source = 'nav';
end
good = find(all(isfinite(fitPosition),2));
inertial = zeros(numel(good),3);
for j = 1:numel(good)
    inertial(j,:) = (rotations(:,:,good(j))'*fitPosition(good(j),:)')';
end
design = [ones(numel(good),1), epochs(good)-tRef];
coefficients = zeros(2,3);
for axis = 1:3
    coefficients(:,axis) = fitL1(design, inertial(:,axis));
end
init.state = [coefficients(1,:)'; coefficients(2,:)'];
init.source = source;
init.snapshotPosition = snapshot;
init.snapshotClockM = snapshotClock;
end

function p = fitL1(A, y)
% 快照解含大量 km 级噪声和偶发发散，用 L1 拟合降低其影响。
p = A\y;
for iteration = 1:100
    w = 1./max(abs(y-A*p),1);
    next = (sqrt(w).*A)\(sqrt(w).*y);
    if max(abs(A*(next-p))) < 1e-3, p = next; break; end
    p = next;
end
end
