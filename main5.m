%% LuGRE 最终算法入口：所有开关及可调参数在 main5_config.m
clear; clc;
projectRoot = fileparts(mfilename('fullpath'));
addpath(projectRoot, fullfile(projectRoot,'function'), fullfile(projectRoot,'reader'));
cfg = main5_config();
results = runLugreAnalysis(cfg);

% 保留常用工作区变量，方便继续绘图；所有钟差统一为米。
NAV = results.nav;
pos = results.dpe.position;
pos2 = results.ls.position;
dt = results.dpe.clockM;
dt2 = results.ls.clockM;
dd = results.dpe.differenceNavM;
dd2 = results.orbit.differenceNavM;
dd3 = results.ls.differenceNavM;
% 真值误差与原有 NAV 差异分开保存，避免混淆参考来源。
if isfield(results,'truth')
    groundTruth = results.truth;
    dpeTruthErrorM = groundTruth.dpe.distanceM;
    lsTruthErrorM = groundTruth.ls.distanceM;
    navTruthErrorM = groundTruth.nav.distanceM;
end
