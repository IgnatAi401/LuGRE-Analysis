function testEpochIntegration()
% TESTEPOCHINTEGRATION 逐历元定位的真实数据与在线 Horizons 集成检查（配置中的默认窗口）。
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root,genpath(fullfile(root,'src')));
originalFolder = pwd;
cleanupFolder = onCleanup(@()cd(originalFolder));
cd(root);
cfg=lugreConfig(); cfg.run.maxEpochs=65; cfg.run.plot=false;
result=runEpochPositioning(cfg);
assert(all(result.dpe.valid) && all(result.ls.valid));
assert(all(isfinite(result.orbit.position),'all'));
assert(size(result.moon,1)==65 && all(abs(result.moon(:,1)-(2444244.5+(result.gpsSeconds-18)/86400))<1e-10));
assert(all(vecnorm(result.moon(:,2:4),2,2)>3e8) && all(vecnorm(result.moon(:,2:4),2,2)<5e8));
fprintf('65-epoch orbit/default-grid integration PASS\n');
% 全部三种求解器组合 × 加权开关；orbit=false 时坏 URL/EOP 不得被调用。
% 真值对比同样需要 EOP 坐标转换，这里关闭，只检查定位流程本身。
for methods=[1 0 1;0 1 1]
    for weighted=[false,true]
        cfg=lugreConfig(); cfg.run.maxEpochs=3; cfg.run.plot=false;
        cfg.enable.dpe=logical(methods(1)); cfg.enable.ls=logical(methods(2));
        cfg.enable.cn0Dpe=weighted; cfg.enable.cn0Ls=weighted;
        cfg.enable.orbit=false; cfg.moon.url='invalid://must-not-fetch';
        cfg.frame.eopSource='must-not-use'; cfg.truth.enabled=false;
        r=runEpochPositioning(cfg);
        assert(isempty(r.moon) && all(isnan(r.orbit.position),'all'));
        assert(all(r.dpe.valid)==cfg.enable.dpe && all(r.ls.valid)==cfg.enable.ls);
        fprintf('switches dpe=%d ls=%d weighted=%d orbit=0 PASS\n',methods,weighted);
    end
end
cfg=lugreConfig(); cfg.run.maxEpochs=3; cfg.run.plot=false;
cfg.enable.dpe=false; cfg.enable.j2=false;
cfg.truth.enabled=false; % 下面只检查定位结果图（3 张），真值图另有 3 张
r=runEpochPositioning(cfg);
assert(all(r.ls.valid) && ~any(r.dpe.valid) && all(isfinite(r.orbit.position),'all'));
fprintf('LS-only orbit with J2 off PASS\n');
originalVisibility = get(groot,'defaultFigureVisible');
cleanupVisibility = onCleanup(@()set(groot,'defaultFigureVisible',originalVisibility));
previousFigures = findall(groot,'Type','figure');
set(groot,'defaultFigureVisible','off'); plotEpochResults(r);
newFigures = setdiff(findall(groot,'Type','figure'),previousFigures);
assert(numel(newFigures)==3); close(newFigures);
fprintf('plot switches PASS\n');
cfg.enable.ls=false;
try
    runEpochPositioning(cfg); error('test:ExpectedFailure','Did not reject disabled solvers');
catch ME
    assert(strcmp(ME.identifier,'LuGRE:NoMethod'));
end
% 静态检查：入口、配置、src 下全部函数与测试。
files=[{'mainEpoch.m','mainBatch.m','lugreConfig.m'}, ...
 cellfun(@(f)erase(f,[root filesep]),listFiles(root),'UniformOutput',false)];
for k=1:numel(files)
    issues=checkcode(files{k},'-id');
    assert(isempty(issues),'LuGRE:CodeCheck','%s 有静态检查诊断。',files{k});
    for j=1:numel(issues)
        fprintf('%s:%d %s %s\n',files{k},issues(j).line,issues(j).id,issues(j).message);
    end
end
fprintf('Integration checks completed.\n');
end

function files = listFiles(root)
entries=[dir(fullfile(root,'src','**','*.m')); dir(fullfile(root,'tests','*.m'))];
files=fullfile({entries.folder},{entries.name});
end
