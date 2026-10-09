function truth = summarizeTruthErrors(names, estimates, position)
% SUMMARIZETRUTHERRORS 各方法相对真值的三维误差，并分解为地心径向（远距离几何的弱方向）与横向。
% names/estimates：方法名与 m×3 ECEF 位置（cell）；position：真值 ECEF。打印汇总表。
truth = struct('position',position);
radialUnit = position./vecnorm(position,2,2);
fprintf('\n%-14s %6s %11s %11s %11s %13s %13s\n','真值对比','历元','中位数/m','RMSE/m','95%/m','径向RMSE/m','横向RMSE/m');
for k = 1:numel(names)
    e = estimates{k} - position;
    radial = sum(e.*radialUnit,2);
    transverse = vecnorm(e - radial.*radialUnit,2,2);
    d = vecnorm(e,2,2);
    good = isfinite(d);
    s = struct('errorM',e,'distanceM',d,'radialM',radial,'transverseM',transverse, ...
        'count',nnz(good),'medianM',median(d(good)),'rmseM',sqrt(mean(d(good).^2)), ...
        'p95M',prctile(d(good),95),'radialRmseM',sqrt(mean(radial(good).^2)), ...
        'transverseRmseM',sqrt(mean(transverse(good).^2)));
    truth.(names{k}) = s;
    fprintf('%-14s %6d %11.0f %11.0f %11.0f %13.0f %13.0f\n',names{k},s.count,s.medianM, ...
        s.rmseM,s.p95M,s.radialRmseM,s.transverseRmseM);
end
end
