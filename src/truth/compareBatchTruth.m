function truth = compareBatchTruth(results)
% COMPAREBATCHTRUTH 批处理各方法、快照 LS 与 NAV 相对真值的误差统计。
% 3D 误差并分解为地心径向（远距离几何的弱方向）与横向。真值只用于事后评估。
cfg = results.config;
[position, valid] = loadTruthEcef(results.gpsSeconds, cfg.truth, cfg.frame);
names = {}; estimates = {};
for method = cfg.batch.methods
    names{end+1} = method{1}; estimates{end+1} = results.(method{1}).positionEcef; %#ok<AGROW>
end
names{end+1} = 'snapshotLs'; estimates{end+1} = results.init.snapshotPosition;
if ~isempty(results.navPosition)
    names{end+1} = 'nav'; estimates{end+1} = results.navPosition;
end
truth = summarizeTruthErrors(names, estimates, position);
truth.valid = valid;
end
