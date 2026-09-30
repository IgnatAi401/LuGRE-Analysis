function testEpoch()
% TESTEPOCH 逐历元方法离线回归：坐标/速度转换、权重、生效开关、无解处理、
% 星历 Toc/Toe、卫星去重、网格精化与钟差消去（需要 Aerospace Toolbox）。
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root,genpath(fullfile(root,'src')));
cfg = lugreConfig();
cfg.dpe.refineStepsM = [];
cfg.dpe.usePredictionCenter = false;
cfg.frame.eopSource = 'manual';
t = 1420996500;
r = [-8e7;4e7;3e7]; v = [100;-200;30];
[ri,vi] = ecefToIcrf(r,v,t,cfg.frame);
[rr,vv] = icrfToEcef(ri,vi,t,cfg.frame);
assert(norm(rr-r)<1e-6 && norm(vv-v)<1e-8);
% 独立检查速度等于位置随时间的变化，避免仅靠正反变换自洽。
h = cfg.frame.derivativeStepS;
rp = icrfToItrfRotation(t+h,cfg.frame)*(ri+h*vi);
rm = icrfToItrfRotation(t-h,cfg.frame)*(ri-h*vi);
assert(norm((rp-rm)/(2*h)-v)<1e-3);
utc = datetime(1980,1,6)+seconds(t-18);
expected = dcmeci2ecef('IAU-2000/2006',datevec(utc),37,0,[0,0]);
assert(norm(icrfToItrfRotation(t,cfg.frame)-expected,'fro')<1e-13);

% 具有足够几何秩的合成观测，不依赖真实数据恰好能收敛。
sat = [2e7 0 0;0 2e7 0;0 0 2e7;-2e7 0 0;0 -2e7 0;0 0 -2e7];
truth = [1000;2000;3000]; bias = 500;
obs.satPosition = sat; obs.range = vecnorm(sat-truth',2,2)+bias;
obs.weight = [0.01;1;1;1;1;1];
s = solveEpochLs(obs,truth+[100;-100;50],bias+100,cfg);
assert(s.valid && norm(s.position-truth)<1e-4 && abs(s.clockM-bias)<1e-4);
noisy = obs; noisy.range(1)=noisy.range(1)+1000;
u = solveEpochLs(noisy,truth,bias,cfg);
cfg.enable.cn0Ls = true;
w = solveEpochLs(noisy,truth,bias,cfg);
assert(u.valid && w.valid && norm(w.position-truth)<norm(u.position-truth));

cfg.dpe.coarseOffsetsM = [-1000,0,1000];
cfg.dpe.fineOffsetsM = cfg.dpe.coarseOffsetsM;
cfg.dpe.clockOffsetsM = [-1000,0,1000]; cfg.dpe.chunkSize = 7;
cfg.enable.orbit = false;
s = solveEpochDpe(obs,truth,bias,nan(3,1),1,cfg);
assert(s.valid && norm(s.position-truth)==0 && s.clockM==bias);
% 同一候选集合在不同分块大小下结果一致。
cfg.dpe.chunkSize = 100;
s2 = solveEpochDpe(obs,truth,bias,nan(3,1),1,cfg);
assert(isequal(s.position,s2.position) && s.cost==s2.cost);
% 轨道先验可被显式开启；关闭时不受强先验参数影响。
cfg.dpe.orbitWeight = 1e6; prior = truth+[1000;0;0];
s = solveEpochDpe(obs,truth,bias,prior,2,cfg);
assert(isequal(s.position,truth));
cfg.enable.orbit = true;
s = solveEpochDpe(obs,truth,bias,prior,2,cfg);
assert(isequal(s.position,prior));

% 直接算两个候选的代价，验证 C/N0 权重确实进入 DPE。
cfg.enable.orbit = false; cfg.enable.cn0Dpe = false;
s1 = solveEpochDpe(noisy,truth,bias,prior,2,cfg);
cfg.enable.cn0Dpe = true;
s2 = solveEpochDpe(noisy,truth,bias,prior,2,cfg);
assert(s1.cost ~= s2.cost);
empty = struct('satPosition',zeros(0,3),'range',[],'weight',[]);
s = solveEpochDpe(empty,truth,bias,prior,1,cfg); assert(~s.valid && all(isnan(s.position)));
s = solveEpochLs(empty,truth,bias,cfg); assert(~s.valid);
% 共用观测模型：Sagnac 开关、CN0 参考尺度、无效伪距筛除。
data = [sat,zeros(6,3),obs.range,[10;30;30;30;30;30]];
cfg.enable.sagnac = false;
off = buildObservations(data,truth,cfg);
assert(isequal(off.satPosition,sat) && abs(off.weight(1)-0.01)<1e-12);
cfg.enable.sagnac = true;
on = buildObservations(data,truth,cfg);
assert(norm(on.satPosition-off.satPosition,'fro')>1);
data(1,7)=NaN; filtered=buildObservations(data,truth,cfg);
assert(numel(filtered.range)==5);

% J2 开关改变加速度但不改变运动学分量。
state = [ri;vi;3.8e8;1e8;1e7;0;1000;0];
cfg.enable.j2 = false; a = orbitDynamics(0,state,t,cfg);
cfg.enable.j2 = true; b = orbitDynamics(0,state,t,cfg);
assert(norm(a(4:6)-b(4:6))>0 && isequal(a(1:3),b(1:3)));
checkMeasurementsAndRefinement(root);
fprintf('testEpoch: all checks passed.\n');
end

function checkMeasurementsAndRefinement(root)
% 卫星去重、缺星历、Toc/Toe 跨周、网格精化与钟差消去。
cfg=lugreConfig(); cfg.enable.orbit=false;
% 相同速度的不同卫星不能合并；同一卫星不同信号只保留一个。
r=struct('rxTime',1,'txTime',0,'satPosition',[2e7;0;0],'satVelocity',zeros(3,1),'prRaw',2e7, ...
 'satClockM',0,'cn0',30,'fdRaw',0,'fdRateRaw',0,'carrierPhase',0,'signalId',0,'svId',1);
raw=repmat(r,1,5); raw(2).signalId=1; raw(3).svId=2;
raw(4).signalId=2; raw(5).cn0=40;
a=selectEpochMeasurements(1,raw); assert(size(a,1)==3 && a(1,8)==40);
raw(1).satPosition(:)=NaN; a=selectEpochMeasurements(1,raw); assert(size(a,1)==3);
warning('off','LuGRE:EphemerisMissing');
p=addSatelliteStates(r,struct([]),struct([]));
warning('on','LuGRE:EphemerisMissing');
assert(all(isnan(p.satVelocity)) && isempty(selectEpochMeasurements(1,p)));
% Toe/Toc 分离：改变钟差参考时间不能改变轨道；跨周参考正确。
[g,~]=readRinexNav(fullfile(root,'data','rinex','BRDC00IGS_R_20250150000_01D_MN.rnx')); e=g(1);
e.M0=0; e.Toe=604790; e.ToeTotalSeconds=604800*2300-10;
e.TocTotalSeconds=e.ToeTotalSeconds+20;
e.ClockDrift=1e-8; e.ClockDriftRate=0;
t=e.ToeTotalSeconds;
[p,v,b]=computeSatellitePvt(e,t,3.986005e14,cfg.constants.omegaEarth,cfg.constants.c);
assert(all(isfinite(v)));
e2=e; e2.TocTotalSeconds=e.TocTotalSeconds+10;
[p2,~,b2]=computeSatellitePvt(e2,t,3.986005e14,cfg.constants.omegaEarth,cfg.constants.c);
assert(norm(p-p2)==0 && abs((b2-b)+10e-8*cfg.constants.c)<1e-7);
h=0.01;
[pp,~,~]=computeSatellitePvt(e,t+h,3.986005e14,cfg.constants.omegaEarth,cfg.constants.c);
[pm,~,~]=computeSatellitePvt(e,t-h,3.986005e14,cfg.constants.omegaEarth,cfg.constants.c);
assert(norm((pp-pm)/(2*h)-v)<0.01);
% 局部精化须降低代价，恢复非粗格点上的合成真值。
sat=[2e7 0 0;0 2e7 0;0 0 2e7;-2e7 0 0;0 -2e7 0;0 0 -2e7];
truth=[120;230;340]; obs.satPosition=sat;
obs.range=vecnorm(sat-truth',2,2)+170; obs.weight=ones(6,1);
cfg.dpe.coarseOffsetsM=-1000:1000:1000; cfg.dpe.clockOffsetsM=-1000:1000:1000;
cfg.dpe.refineStepsM=[]; coarse=solveEpochDpe(obs,zeros(3,1),0,nan(3,1),1,cfg);
checkClockElimination(obs,cfg);
cfg.dpe.refineStepsM=[100 10]; fine=solveEpochDpe(obs,zeros(3,1),0,nan(3,1),1,cfg);
assert(fine.valid && fine.cost<coarse.cost && norm(fine.position-truth)<1e-6);
cfg.enable.orbit=true; cfg.dpe.usePredictionCenter=true;
s=solveEpochDpe(obs,[1e6;0;0],170,truth,2,cfg);
assert(s.valid && norm(s.position-truth)<1e-6);
cfg.enable.orbit=false;
obs.satPosition=repmat(sat(1,:),6,1);
s=solveEpochDpe(obs,zeros(3,1),0,nan(3,1),1,cfg); assert(~s.valid);
end

function checkClockElimination(obs,cfg)
% 非均匀钟差网格、等权/加权、钟差越界、分块和平局均与穷举比较。
cfg.dpe.clockOffsetsM=[-900 -250 0 130 700];
cfg.dpe.chunkSize=7;
for weighted=[false true]
    cfg.enable.cn0Dpe=weighted;
    obs.weight=[0.2;0.5;1;2;3;4];
    for shift=[-5000 0 5000]
        shifted=obs; shifted.range=obs.range+shift;
        s=solveEpochDpe(shifted,zeros(3,1),0,nan(3,1),1,cfg);
        offsets=cfg.dpe.coarseOffsetsM;
        [x,y,z,b]=ndgrid(offsets,offsets,offsets,cfg.dpe.clockOffsetsM);
        pos=[x(:) y(:) z(:)]; cost=zeros(size(x(:)));
        w=ones(6,1); if weighted, w=obs.weight; end
        for j=1:6
            cost=cost+w(j)*abs(shifted.range(j)-vecnorm(pos-obs.satPosition(j,:),2,2)-b(:));
        end
        [best,k]=min(cost);
        assert(abs(s.cost-best)<1e-7 && isequal(s.position,pos(k,:)') && s.clockM==b(k));
    end
end
end
