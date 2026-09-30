%% LuGRE 多历元批处理 DPE 入口：参数在 main5_config.m 的 cfg.data 与 cfg.batch
clear; clc;
projectRoot = fileparts(mfilename('fullpath'));
addpath(projectRoot, fullfile(projectRoot,'function'), fullfile(projectRoot,'reader'));
cfg = main5_config();
batch = runLugreBatch(cfg);
