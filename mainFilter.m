%% LuGRE 扩展卡尔曼滤波入口：ekfOrbit（地月动力学 EKF）/ ekfCvBlend（匀速 EKF + 轨道后加权），参数在 lugreConfig.m 的 cfg.filter
clear; clc;
projectRoot = fileparts(mfilename('fullpath'));
addpath(projectRoot, genpath(fullfile(projectRoot,'src')));
cfg = lugreConfig();
ekf = runFilterPositioning(cfg);
