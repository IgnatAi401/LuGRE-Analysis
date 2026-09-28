function results = compareLugreTruth(results, options)
% COMPARELUGRETRUTH 按 GPST 对齐地心 J2000 km 真值，转换至 ITRF m。
% 真值只用于事后评估，不参与定位或轨道反馈。
header = string(readcell(options.file,'Range','C1:H1'));
assert(isequal(header([1,4:6]), ["gps_seconds", ...
    "lugre_rel_earth_x_j2000_km","lugre_rel_earth_y_j2000_km", ...
    "lugre_rel_earth_z_j2000_km"]),'LuGRE:TruthSchema','真值表列名或单位不符。');
data = readmatrix(options.file,'Range','C2:H1048576');
data = data(all(isfinite(data(:,[1,4:6])),2),[1,4:6]);
query = results.gpsSeconds(:);
data = data(data(:,1)>=min(query)-options.maxGapS & ...
    data(:,1)<=max(query)+options.maxGapS,:);
assert(~isempty(data),'LuGRE:TruthCoverage','真值没有覆盖结果时间。');
data = sortrows(data,1);
[time, first, group] = unique(data(:,1));
assert(isequal(data(:,2:4),data(first(group),2:4)), ...
    'LuGRE:TruthDuplicate','同一真值时间存在不同坐标。');
position = data(first,2:4)*1000;
assert(numel(time)>=2,'LuGRE:TruthCoverage','至少需要两个真值时间点。');
inertial = interp1(time,position,query,'linear',NaN);
left = interp1(time,(1:numel(time))',query,'previous',NaN);
right = interp1(time,(1:numel(time))',query,'next',NaN);
valid = isfinite(left) & isfinite(right);
idx = find(valid);
valid(idx) = time(right(idx))-time(left(idx)) <= options.maxGapS;
assert(any(valid),'LuGRE:TruthCoverage','结果时间内没有可用的连续真值。');
ecef = nan(numel(query),3);
for k = find(valid)'
    rotation = lugreFrameRotation(query(k),results.config.frame);
    ecef(k,:) = (rotation*inertial(k,:)')';
end
results.truth = struct('source',options.file,'position',ecef,'valid',valid, ...
    'gpsSeconds',query,'alignment','linear GPST; no extrapolation or gap bridging');
names = {'dpe','ls','orbit','nav'};
for k = 1:numel(names)
    name = names{k};
    if strcmp(name,'nav')
        p = results.referencePosition;
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
