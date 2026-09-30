%% LuGRE 逐历元定位入口（最终版算法）：网格 DPE + LS 对照，参数在 lugreConfig.m
clear; clc;
projectRoot = fileparts(mfilename('fullpath'));
addpath(projectRoot, genpath(fullfile(projectRoot,'src')));
cfg = lugreConfig();
results = runEpochPositioning(cfg);

% 常用工作区变量，方便继续绘图；钟差统一为米。
dpePosition = results.dpe.position;
lsPosition = results.ls.position;
orbitPosition = results.orbit.position;
dpeClockM = results.dpe.clockM;
lsClockM = results.ls.clockM;
if isfield(results,'truth')
    dpeTruthErrorM = results.truth.dpe.distanceM;
    lsTruthErrorM = results.truth.ls.distanceM;
    navTruthErrorM = results.truth.nav.distanceM;
end
