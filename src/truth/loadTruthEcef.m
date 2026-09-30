function [ecef, valid] = loadTruthEcef(gpsSeconds, options, frame)
% LOADTRUTHECEF 按 GPST 线性插值地心 J2000 km 真值，转换至 ITRF m。
% 不外推，不跨越超过 options.maxGapS 的缺口；真值只用于事后评估。
data = readTruthTable(options);
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
rotations = icrfToItrfRotation(query(valid),frame);
ecef(valid,:) = squeeze(pagemtimes(rotations,permute(inertial(valid,:),[2 3 1])))';
end

function data = readTruthTable(options)
% 读取 [GPS 秒, J2000 x y z (km)]。xlsx 解析每次约 5 s，结果缓存到 options.cacheFile（.mat）；
% 源文件大小或修改时间变化时自动重新解析。未配置 cacheFile 时不缓存。
source = dir(options.file);
assert(isscalar(source),'LuGRE:TruthMissing','找不到真值文件 %s。',options.file);
stamp = [source.bytes, source.datenum];
useCache = isfield(options,'cacheFile') && ~isempty(options.cacheFile);
if useCache && isfile(options.cacheFile)
    cached = load(options.cacheFile,'stamp','data');
    if isequal(cached.stamp,stamp), data = cached.data; return; end
end
table = readtable(options.file,'Range','C:H','VariableNamingRule','preserve');
assert(isequal(string(table.Properties.VariableNames([1,4:6])), ["gps_seconds", ...
    "lugre_rel_earth_x_j2000_km","lugre_rel_earth_y_j2000_km", ...
    "lugre_rel_earth_z_j2000_km"]),'LuGRE:TruthSchema','真值表列名或单位不符。');
data = table{:,[1,4:6]};
data = data(all(isfinite(data),2),:);
if useCache
    folder = fileparts(options.cacheFile);
    if ~isfolder(folder), mkdir(folder); end
    save(options.cacheFile,'stamp','data');
end
end
