function [summary, status, wide] = summarizeMethodCompare(folder)
% SUMMARIZEMETHODCOMPARE 汇总 compareMethods 每个 OP 的 .mat：
% all_summary.csv（各 OP 各方法的统计，长表）、all_status.csv（OP × 方法的运行状态）、
% all_median.csv（OP × 方法：自身历元三维误差中位数，km；另附共同历元中位数表 all_commonMedian.csv）。
% 早期运行的 statusTable 对未收敛的批处理也记为 ok，这里按 info.converged 统一改为 ok (not converged)。
% 各 OP 的方法集合可能不同（如后加的 ekfLunarStatic 只在部分 .mat 中），取并集，缺的状态记为 not run。
root = fileparts(fileparts(mfilename('fullpath')));
if nargin < 1, folder = fullfile(root,'results','method_compare_all'); end
files = dir(fullfile(folder,'OP*_*.mat'));
files = files(~contains({files.name},'_epochL1'));
% 按 OP 编号排序
key = cellfun(@(n)sscanf(n,'OP%d_%d')'*[100;1], {files.name});
[~,order] = sort(key); files = files(order);
loaded = arrayfun(@(f)load(fullfile(folder,f.name),'info','stats','statusTable'), files);
% 方法顺序取方法最多的那个 OP，其余补在后面。
lists = arrayfun(@(s){s.statusTable.method'}, loaded);
[~,longest] = max(cellfun(@numel, lists));
allMethods = unique([lists{longest}, lists{:}], 'stable');
summary = table(); status = table(); median3d = table(); commonMedian = table();
for k = 1:numel(files)
    s = loaded(k);
    info = s.info;
    s.statusTable = markNotConverged(s.statusTable, info);
    missing = setdiff(allMethods, s.statusTable.method', 'stable');
    if ~isempty(missing)
        s.statusTable = [s.statusTable; cell2table([missing', repmat({'not run'},numel(missing),1)], ...
            'VariableNames',{'method','status'})];
    end
    [~,order] = ismember(allMethods, s.statusTable.method');
    s.statusTable = s.statusTable(order,:);
    summary = [summary; s.stats]; %#ok<AGROW>
    meta = table(string(info.window), info.rawEpochs, info.epochsWith4, info.navEpochs, info.meanRangeRe, ...
        'VariableNames',{'window','rawEpochs','epochsWith4','navEpochs','meanRangeRe'});
    st = cell2table(s.statusTable.status','VariableNames',s.statusTable.method');
    status = [status; [meta st]]; %#ok<AGROW>
    methods = s.statusTable.method';
    own = nan(1,numel(methods)); both = own;
    for j = 1:numel(methods)
        hit = strcmp(s.stats.method, methods{j});
        if any(hit), own(j) = s.stats.medianM(hit)/1000; both(j) = s.stats.commonMedianM(hit)/1000; end
    end
    median3d = [median3d; [meta array2table(own,'VariableNames',methods)]]; %#ok<AGROW>
    commonMedian = [commonMedian; [meta array2table(both,'VariableNames',methods)]]; %#ok<AGROW>
end
writetable(summary, fullfile(folder,'all_summary.csv'));
writetable(status, fullfile(folder,'all_status.csv'));
writetable(median3d, fullfile(folder,'all_median.csv'));
writetable(commonMedian, fullfile(folder,'all_commonMedian.csv'));
wide = median3d;
disp(status(:,1:5)); disp(median3d);
end

function statusTable = markNotConverged(statusTable, info)
if ~isfield(info,'converged'), return; end
map = struct('dpeGrid','batchDpeGrid','dpeIrls','batchDpeIrls','ls','batchLs');
for name = fieldnames(map)'
    row = strcmp(statusTable.method, map.(name{1}));
    if any(row) && strcmp(statusTable.status{row},'ok') && ~info.converged.(name{1})
        statusTable.status{row} = 'ok (not converged)';
    end
end
end
