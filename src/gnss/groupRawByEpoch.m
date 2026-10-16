function groups = groupRawByEpoch(raw, epochs)
% GROUPRAWBYEPOCH 按接收历元把 RAW 记录分组，groups{k} 为 rxTime 在 epochs(k)±1e-6 s 内的记录（保持原顺序）。
% 供 selectEpochMeasurements(epochs(k), groups{k}) 使用，与直接传入全部 RAW 结果相同，
% 但避免每个历元都扫描全部记录（长窗口上耗时随历元数平方增长）。epochs 须严格递增且间隔 > 2e-6 s。
epochs = epochs(:);
groups = cell(numel(epochs),1);
if isempty(raw) || isempty(epochs), return; end
assert(all(diff(epochs) > 2e-6), 'LuGRE:EpochOrder', '历元必须严格递增。');
edges = reshape([epochs-1e-6, epochs+1e-6]', [], 1);
[~,~,bin] = histcounts([raw.rxTime]', edges);
inside = mod(bin,2) == 1; % 奇数 bin 落在某个历元的容差区间内
groups(:) = {raw([])};
index = find(inside);
members = accumarray((bin(index)+1)/2, index, [numel(epochs),1], @(v){sort(v)}, {zeros(0,1)});
for k = 1:numel(epochs)
    if ~isempty(members{k}), groups{k} = raw(members{k}); end
end
end
