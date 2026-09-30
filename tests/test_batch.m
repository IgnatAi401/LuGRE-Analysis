function test_batch()
% TEST_BATCH 批处理 DPE/LS 的合成数据检查（离线，不需要 Aerospace Toolbox）。
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'function'),fullfile(root,'reader'));
cfg = main5_config();
cfg.enable.j2 = false;
cfg.orbit.relativeTolerance = 1e-11;
cfg.batch.clockModel = 'polynomial'; cfg.batch.clockDegree = 2;
cfg.batch.signalBias = true; % 合成场景每种信号有 3 颗星，偏差可估；真实数据默认关闭

% 轨道积分：正向再反向回到原点应一致。
tRef = 1421300000;
moon0 = [3.8e8;0;0;0;1000;0];
state = [8e7;5e7;2e7;-800;1200;300];
[r1,v1,moon1] = lugrePropagateArc(state,moon0,tRef,tRef+[-300;0;600],cfg);
assert(isequal(r1(2,:)',state(1:3)) && all(isfinite(r1),'all'));
back = lugrePropagateArc([r1(3,:)';v1(3,:)'],moon1(3,:)',tRef+600,tRef,cfg);
assert(norm(back-state(1:3)') < 1, '正反积分不一致：%.3f m', norm(back-state(1:3)'));

% 合成场景：6 颗 GNSS 星，每历元随机 1~5 颗可用，钟差二次变化，两类信号有常值偏差。
rng(7);
epochs = tRef + (-600:600)';
m = numel(epochs);
truthIcrf = lugrePropagateArc(state,moon0,tRef,epochs,cfg);
rotations = repmat(eye(3),1,1,m);
radius = 2.656e7; n = 1.4585e-4;
planes = [0 55 110 165 220 275]*pi/180; phase = (0:5)*1.1;
clock = @(t) -1e5 - 60*(t-tRef) + 0.002*(t-tRef).^2;
rows = {};
for k = 1:m
    visible = randperm(6, randi(5));
    for j = visible
        angle = phase(j) + n*(epochs(k)-tRef);
        node = [cos(planes(j)); sin(planes(j)); 0];
        normal = [0;0;1]; inPlane = cross(normal,node)*cosd(55) + normal*sind(55);
        sat = radius*(cos(angle)*node + sin(angle)*inPlane);
        signal = 0; if j > 3, signal = 2; end
        range = norm(truthIcrf(k,:)'-sat) + clock(epochs(k)) + 7*(signal==2) + 3*randn;
        rows{end+1} = [k, sat', range, signal]; %#ok<AGROW>
    end
end
data = cell2mat(rows');
obs = struct('epochIndex',data(:,1),'satPosition',data(:,2:4),'range',data(:,5), ...
    'weight',ones(size(data,1),1),'signalId',data(:,6));
start = state + [5e3;-4e3;3e3;5;-5;3];
cfg.batch.gridHalfWidthM = 20000; cfg.batch.gridStepM = 500;

clean = struct();
for method = {'dpe','ls'}
    s = lugreBatchSolve(obs,epochs,rotations,start,moon0,tRef,method{1},cfg);
    err = vecnorm(s.positionEcef-truthIcrf,2,2);
    assert(s.converged && max(err) < 5*max(s.formalStd3dM) + 1, ...
        '%s 无粗差时误差 %.1f m 超过形式 5σ %.1f m', method{1}, max(err), 5*max(s.formalStd3dM));
    % 参考信号取观测最多者；signal 2 相对 signal 0 的真值偏差为 +7 m。
    expected = 7*sign(s.biasSignals - s.referenceSignal);
    assert(isscalar(s.biasSignals) && abs(s.signalBiasM - expected) < 3, ...
        '信号偏差 %.2f m（期望 %.0f m）', s.signalBiasM, expected);
    assert(max(abs(s.clockM - clock(epochs))) < 5*max(s.formalStd3dM) + 5);
    clean.(method{1}) = max(err);
end

% 2% 观测加 +3 km 粗差：L1 代价的 DPE 应基本不受影响，且优于 L2 的 LS。
bad = obs;
hit = rand(size(bad.range)) < 0.02;
bad.range(hit) = bad.range(hit) + 3000;
dpe = lugreBatchSolve(bad,epochs,rotations,start,moon0,tRef,'dpe',cfg);
ls = lugreBatchSolve(bad,epochs,rotations,start,moon0,tRef,'ls',cfg);
dpeErr = max(vecnorm(dpe.positionEcef-truthIcrf,2,2));
lsErr = max(vecnorm(ls.positionEcef-truthIcrf,2,2));
assert(dpeErr < 2*clean.dpe + 5 && dpeErr < lsErr, ...
    '粗差下 DPE %.1f m、LS %.1f m（无粗差 DPE %.1f m）', dpeErr, lsErr, clean.dpe);
fprintf('test_batch: 无粗差 DPE %.1f m / LS %.1f m；粗差下 DPE %.1f m / LS %.1f m\n', ...
    clean.dpe, clean.ls, dpeErr, lsErr);

% 样条钟差：二次钟差也能被表示，且样条基函数构成单位分解（常值钟差可精确表示）。
splineCfg = cfg; splineCfg.batch.clockModel = 'spline'; splineCfg.batch.clockKnotSpacingS = 300;
s = lugreBatchSolve(obs,epochs,rotations,start,moon0,tRef,'ls',splineCfg);
assert(s.converged && max(vecnorm(s.positionEcef-truthIcrf,2,2)) < 5*max(s.formalStd3dM) + 1);
assert(numel(s.clockCoefficientsM) == ceil(1200/300) + 3);

% 关闭信号偏差估计时不能生成全零列。
cfg.batch.signalBias = false;
s = lugreBatchSolve(obs,epochs,rotations,start,moon0,tRef,'ls',cfg);
assert(isempty(s.biasSignals) && s.parameterCount == 6 + cfg.batch.clockDegree + 1);

% 自动选星历：日期匹配，缺文件明确报错。
cfg.data.rinex = 'auto'; cfg.data.rinexDir = root;
files = selectRinexFiles(cfg, 1420991587);
assert(endsWith(files{1},'BRDC00IGS_R_20250150000_01D_MN.rnx'));
try
    selectRinexFiles(cfg, 1421117218); % 2025-01-17，本项目没有该日星历
    error('test:NoError','应报缺少星历');
catch ME
    assert(strcmp(ME.identifier,'LuGRE:RinexMissing'));
end
% 星历龄期检查：10 天外的星历必须被拒绝。
[g,~] = read_gnss_rinex(files{1});
r = struct('rxTime',0,'txTime',g(1).ToeTotalSeconds+10*86400,'signalId',0,'svId',g(1).PRN);
warning('off','LuGRE:EphemerisStale'); cleanup = onCleanup(@()warning('on','LuGRE:EphemerisStale'));
stale = getPosRAW(r,g,struct([]));
fresh = getPosRAW(r,g,struct([]),Inf);
assert(all(isnan(stale.ECEF)) && all(isfinite(fresh.ECEF)));
fprintf('test_batch: all checks passed.\n');
end
