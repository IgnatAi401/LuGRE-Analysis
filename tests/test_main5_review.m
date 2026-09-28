function test_main5_review()
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'reader'),fullfile(root,'function'));
cfg=main5_config(); cfg.enable.orbit=false;
% 相同速度的不同卫星不能合并；同一卫星不同信号只保留一个。
r=struct('rxTime',1,'txTime',0,'ECEF',[2e7;0;0],'Vel',zeros(3,1),'prRaw',2e7, ...
 'clockBias',0,'cn0',30,'fdRaw',0,'fdRateRaw',0,'carrierPhase',0,'signalId',0,'svId',1);
raw=repmat(r,1,5); raw(2).signalId=1; raw(3).svId=2;
raw(4).signalId=2; raw(5).cn0=40;
a=getPos(1,raw); assert(size(a,1)==3 && a(1,8)==40);
raw(1).ECEF(:)=NaN; a=getPos(1,raw); assert(size(a,1)==3);
p=getPosRAW(r,struct([]),struct([]));
assert(all(isnan(p.Vel)) && isempty(getPos(1,p)));
% Toe/Toc 分离：改变钟差参考时间不能改变轨道；跨周参考正确。
[g,~]=read_gnss_rinex(cfg.data.rinex); e=g(1);
e.M0=0; e.Toe=604790; e.ToeTotalSeconds=604800*2300-10;
e.TocTotalSeconds=e.ToeTotalSeconds+20; e.GPSTotalSeconds=e.TocTotalSeconds;
e.ClockDrift=1e-8; e.ClockDriftRate=0;
t=e.ToeTotalSeconds;
[p,v,b]=calculateSatPVT(e,t,3.986005e14,cfg.constants.omegaEarth,cfg.constants.c);
assert(all(isfinite(v)));
e2=e; e2.TocTotalSeconds=e.TocTotalSeconds+10;
[p2,~,b2]=calculateSatPVT(e2,t,3.986005e14,cfg.constants.omegaEarth,cfg.constants.c);
assert(norm(p-p2)==0 && abs((b2-b)+10e-8*cfg.constants.c)<1e-7);
h=0.01;
[pp,~,~]=calculateSatPVT(e,t+h,3.986005e14,cfg.constants.omegaEarth,cfg.constants.c);
[pm,~,~]=calculateSatPVT(e,t-h,3.986005e14,cfg.constants.omegaEarth,cfg.constants.c);
assert(norm((pp-pm)/(2*h)-v)<0.01);
legacy=rmfield(e,{'ToeTotalSeconds','TocTotalSeconds'});
[pl,~,~]=calculateSatPVT(legacy,t,3.986005e14,cfg.constants.omegaEarth,cfg.constants.c);
assert(norm(pl-p)<1e-6);
% 局部精化须降低代价，恢复非粗格点上的合成真值。
sat=[2e7 0 0;0 2e7 0;0 0 2e7;-2e7 0 0;0 -2e7 0;0 0 -2e7];
truth=[120;230;340]; obs.position=sat;
obs.range=vecnorm(sat-truth',2,2)+170; obs.weight=ones(6,1);
cfg.dpe.coarseOffsetsM=-1000:1000:1000; cfg.dpe.clockOffsetsM=-1000:1000:1000;
cfg.dpe.refineStepsM=[]; coarse=lugreDpe(obs,zeros(3,1),0,nan(3,1),1,cfg);
checkClockElimination(obs,cfg);
cfg.dpe.refineStepsM=[100 10]; fine=lugreDpe(obs,zeros(3,1),0,nan(3,1),1,cfg);
assert(fine.valid && fine.cost<coarse.cost && norm(fine.position-truth)<1e-6);
cfg.enable.orbit=true; cfg.dpe.usePredictionCenter=true;
s=lugreDpe(obs,[1e6;0;0],170,truth,2,cfg);
assert(s.valid && norm(s.position-truth)<1e-6);
cfg.enable.orbit=false;
obs.position=repmat(sat(1,:),6,1);
s=lugreDpe(obs,zeros(3,1),0,nan(3,1),1,cfg); assert(~s.valid);
fprintf('test_main5_review: all checks passed.\n');
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
        s=lugreDpe(shifted,zeros(3,1),0,nan(3,1),1,cfg);
        offsets=cfg.dpe.coarseOffsetsM;
        [x,y,z,b]=ndgrid(offsets,offsets,offsets,cfg.dpe.clockOffsetsM);
        pos=[x(:) y(:) z(:)]; cost=zeros(size(x(:)));
        w=ones(6,1); if weighted, w=obs.weight; end
        for j=1:6
            cost=cost+w(j)*abs(shifted.range(j)-vecnorm(pos-obs.position(j,:),2,2)-b(:));
        end
        [best,k]=min(cost);
        assert(abs(s.cost-best)<1e-7 && isequal(s.position,pos(k,:)') && s.clockM==b(k));
    end
end
end
