function test_main5_integration()
% TEST_MAIN5_INTEGRATION 真实 OP1_0 数据与在线 Horizons 集成检查。
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'function'));
originalFolder = pwd;
cleanupFolder = onCleanup(@()cd(originalFolder)); %#ok<NASGU>
cd(root);
cfg=main5_config(); cfg.run.maxEpochs=65; cfg.run.plot=false;
result=runLugreAnalysis(cfg);
assert(all(result.dpe.valid) && all(result.ls.valid));
assert(all(isfinite(result.orbit.position),'all'));
assert(size(result.moon,1)==65 && all(abs(result.moon(:,1)-(2444244.5+(result.gpsSeconds-18)/86400))<1e-10));
assert(all(vecnorm(result.moon(:,2:4),2,2)>3e8) && all(vecnorm(result.moon(:,2:4),2,2)<5e8));
fprintf('65-epoch orbit/default-grid integration PASS\n');
% 全部三种求解器组合 × 加权开关；orbit=false 时坏 URL/EOP 不得被调用。
for methods=[1 0 1;0 1 1]
    for weighted=[false,true]
        cfg=main5_config(); cfg.run.maxEpochs=3; cfg.run.plot=false;
        cfg.enable.dpe=logical(methods(1)); cfg.enable.ls=logical(methods(2));
        cfg.enable.cn0Dpe=weighted; cfg.enable.cn0Ls=weighted;
        cfg.enable.orbit=false; cfg.moon.url='invalid://must-not-fetch';
        cfg.frame.eopSource='must-not-use';
        r=runLugreAnalysis(cfg);
        assert(isempty(r.moon) && all(isnan(r.orbit.position),'all'));
        assert(all(r.dpe.valid)==cfg.enable.dpe && all(r.ls.valid)==cfg.enable.ls);
        fprintf('switches dpe=%d ls=%d weighted=%d orbit=0 PASS\n',methods,weighted);
    end
end
cfg=main5_config(); cfg.run.maxEpochs=3; cfg.run.plot=false;
cfg.enable.dpe=false; cfg.enable.j2=false;
r=runLugreAnalysis(cfg);
assert(all(r.ls.valid) && ~any(r.dpe.valid) && all(isfinite(r.orbit.position),'all'));
fprintf('LS-only orbit with J2 off PASS\n');
originalVisibility = get(groot,'defaultFigureVisible');
cleanupVisibility = onCleanup(@()set(groot,'defaultFigureVisible',originalVisibility)); %#ok<NASGU>
previousFigures = findall(groot,'Type','figure');
set(groot,'defaultFigureVisible','off'); plotLugreResults(r);
newFigures = setdiff(findall(groot,'Type','figure'),previousFigures);
assert(numel(newFigures)==3); close(newFigures);
fprintf('plot switches PASS\n');
cfg.enable.ls=false;
try
    runLugreAnalysis(cfg); error('test:ExpectedFailure','Did not reject disabled solvers');
catch ME
    assert(strcmp(ME.identifier,'LuGRE:NoMethod'));
end
files={'main5.m','main5_config.m','runLugreAnalysis.m','function/lugreFrameRotation.m', ...
 'function/getMoonEph.m','function/ecef_to_j2000.m','function/j2000_to_ecef.m', ...
 'function/lugreDpe.m','function/lugreLs.m','function/lugreDynamics.m', ...
 'function/lugreObservations.m','function/plotLugreResults.m','tests/test_main5.m'};
for k=1:numel(files)
    issues=checkcode(files{k},'-id');
    assert(isempty(issues),'LuGRE:CodeCheck','%s 有静态检查诊断。',files{k});
    for j=1:numel(issues)
        fprintf('%s:%d %s %s\n',files{k},issues(j).line,issues(j).id,issues(j).message);
    end
end
fprintf('Integration checks completed.\n');
end
