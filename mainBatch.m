%% LuGRE 多历元批处理入口：dpeIrls / dpeGrid 与同模型 LS 对照，参数在 lugreConfig.m 的 cfg.batch
clear; clc;
projectRoot = fileparts(mfilename('fullpath'));
addpath(projectRoot, genpath(fullfile(projectRoot,'src')));
cfg = lugreConfig();
batch = runBatchPositioning(cfg);
