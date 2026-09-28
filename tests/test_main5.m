function test_main5()
% TEST_MAIN5 离线回归：坐标/速度转换、权重、生效开关、无解处理。
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'function'));
cfg = main5_config();
cfg.dpe.refineStepsM = [];
cfg.dpe.usePredictionCenter = false;
cfg.frame.eopSource = 'manual';
t = 1420996500;
r = [-8e7;4e7;3e7]; v = [100;-200;30];
[ri,vi] = ecef_to_j2000(r,v,t,cfg.frame);
[rr,vv] = j2000_to_ecef(ri,vi,t,cfg.frame);
assert(norm(rr-r)<1e-6 && norm(vv-v)<1e-8);
% 独立检查速度等于位置随时间的变化，避免仅靠正反变换自洽。
h = cfg.frame.derivativeStepS;
rp = lugreFrameRotation(t+h,cfg.frame)*(ri+h*vi);
rm = lugreFrameRotation(t-h,cfg.frame)*(ri-h*vi);
assert(norm((rp-rm)/(2*h)-v)<1e-3);
utc = datetime(1980,1,6)+seconds(t-18);
expected = dcmeci2ecef('IAU-2000/2006',datevec(utc),37,0,[0,0]);
assert(norm(lugreFrameRotation(t,cfg.frame)-expected,'fro')<1e-13);

% 具有足够几何秩的合成观测，不依赖真实数据恰好能收敛。
sat = [2e7 0 0;0 2e7 0;0 0 2e7;-2e7 0 0;0 -2e7 0;0 0 -2e7];
truth = [1000;2000;3000]; bias = 500;
obs.position = sat; obs.range = vecnorm(sat-truth',2,2)+bias;
obs.weight = [0.01;1;1;1;1;1];
s = lugreLs(obs,truth+[100;-100;50],bias+100,cfg);
assert(s.valid && norm(s.position-truth)<1e-4 && abs(s.clockM-bias)<1e-4);
noisy = obs; noisy.range(1)=noisy.range(1)+1000;
u = lugreLs(noisy,truth,bias,cfg);
cfg.enable.cn0Ls = true;
w = lugreLs(noisy,truth,bias,cfg);
assert(u.valid && w.valid && norm(w.position-truth)<norm(u.position-truth));

cfg.dpe.coarseOffsetsM = [-1000,0,1000];
cfg.dpe.fineOffsetsM = cfg.dpe.coarseOffsetsM;
cfg.dpe.clockOffsetsM = [-1000,0,1000]; cfg.dpe.chunkSize = 7;
cfg.enable.orbit = false;
s = lugreDpe(obs,truth,bias,nan(3,1),1,cfg);
assert(s.valid && norm(s.position-truth)==0 && s.clockM==bias);
% 同一候选集合在不同分块大小下结果一致。
cfg.dpe.chunkSize = 100;
s2 = lugreDpe(obs,truth,bias,nan(3,1),1,cfg);
assert(isequal(s.position,s2.position) && s.cost==s2.cost);
% 轨道先验可被显式开启；关闭时不受强先验参数影响。
cfg.dpe.orbitWeight = 1e6; prior = truth+[1000;0;0];
s = lugreDpe(obs,truth,bias,prior,2,cfg);
assert(isequal(s.position,truth));
cfg.enable.orbit = true;
s = lugreDpe(obs,truth,bias,prior,2,cfg);
assert(isequal(s.position,prior));

% 直接算两个候选的代价，验证 C/N0 权重确实进入 DPE。
cfg.enable.orbit = false; cfg.enable.cn0Dpe = false;
s1 = lugreDpe(noisy,truth,bias,prior,2,cfg);
cfg.enable.cn0Dpe = true;
s2 = lugreDpe(noisy,truth,bias,prior,2,cfg);
assert(s1.cost ~= s2.cost);
empty = struct('position',zeros(0,3),'range',[],'weight',[]);
s = lugreDpe(empty,truth,bias,prior,1,cfg); assert(~s.valid && all(isnan(s.position)));
s = lugreLs(empty,truth,bias,cfg); assert(~s.valid);
% 共用观测模型：Sagnac 开关、CN0 参考尺度、无效伪距筛除。
data = [sat,zeros(6,3),obs.range,[10;30;30;30;30;30]];
cfg.enable.sagnac = false;
off = lugreObservations(data,truth,cfg);
assert(isequal(off.position,sat) && abs(off.weight(1)-0.01)<1e-12);
cfg.enable.sagnac = true;
on = lugreObservations(data,truth,cfg);
assert(norm(on.position-off.position,'fro')>1);
data(1,7)=NaN; filtered=lugreObservations(data,truth,cfg);
assert(numel(filtered.range)==5);

% J2 开关改变加速度但不改变运动学分量。
state = [ri;vi;3.8e8;1e8;1e7;0;1000;0];
cfg.enable.j2 = false; a = lugreDynamics(0,state,t,cfg);
cfg.enable.j2 = true; b = lugreDynamics(0,state,t,cfg);
assert(norm(a(4:6)-b(4:6))>0 && isequal(a(1:3),b(1:3)));
test_main5_review();
fprintf('test_main5: all checks passed.\n');
end
