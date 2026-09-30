function results = compareEpochTruth(results, options)
% COMPAREEPOCHTRUTH 逐历元 DPE、LS、轨道预测与 NAV 相对真值的误差统计。
% 真值只用于事后评估，不参与定位或轨道反馈。
query = results.gpsSeconds(:);
[ecef, valid] = loadTruthEcef(query, options, results.config.frame);
results.truth = struct('source',options.file,'position',ecef,'valid',valid, ...
    'gpsSeconds',query,'alignment','linear GPST; no extrapolation or gap bridging');
names = {'dpe','ls','orbit','nav'};
for k = 1:numel(names)
    name = names{k};
    if strcmp(name,'nav')
        p = results.navPosition;
    else
        p = results.(name).position;
        if isfield(results.(name),'valid'), p(~results.(name).valid,:) = NaN; end
    end
    errorXYZ = p-ecef;
    distance = vecnorm(errorXYZ,2,2);
    good = isfinite(distance);
    results.truth.(name) = struct('errorM',errorXYZ,'distanceM',distance, ...
        'count',nnz(good),'rmseM',sqrt(mean(distance(good).^2)), ...
        'medianM',median(distance(good)));
    fprintf('%s truth: %d/%d, 3D RMSE %.3f m, median %.3f m\n', ...
        upper(name),nnz(good),numel(query),results.truth.(name).rmseM,results.truth.(name).medianM);
end
if any(~valid)
    warning('LuGRE:TruthCoverage','%d 个历元缺少连续真值，误差保留 NaN。',nnz(~valid));
end
end
