function [ecef, valid] = loadTruthEcef(gpsSeconds, options, frame)
% LOADTRUTHECEF 按 GPST 线性插值地心 J2000 km 真值，转换至 ITRF m。
% 不外推，不跨越超过 options.maxGapS 的缺口；真值只用于事后评估。
header = string(readcell(options.file,'Range','C1:H1'));
assert(isequal(header([1,4:6]), ["gps_seconds", ...
    "lugre_rel_earth_x_j2000_km","lugre_rel_earth_y_j2000_km", ...
    "lugre_rel_earth_z_j2000_km"]),'LuGRE:TruthSchema','真值表列名或单位不符。');
data = readmatrix(options.file,'Range','C2:H1048576');
data = data(all(isfinite(data(:,[1,4:6])),2),[1,4:6]);
query = gpsSeconds(:);
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
    rotation = icrfToItrfRotation(query(k),frame);
    ecef(k,:) = (rotation*inertial(k,:)')';
end
end
