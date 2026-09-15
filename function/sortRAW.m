function [RAW_sorted] = sortRAW(NAV, RAW)

%   1. 判断RAW中每组数据的rxTime是否在NAV的rxTime中出现
%   2. 未出现的组被剔除，出现的组保留
%   3. 为保留的组计算 txTime = rxTime - prRaw/光速

    % 光速常量 (m/s)
    c = 299792458;
    
    % 提取NAV中所有的rxTime值
    NAV_rxTime = [NAV.rxTime];
    
    % 提取RAW中所有的rxTime值
    RAW_rxTime = [RAW.rxTime];
    
    % 找出RAW中rxTime在NAV中出现过的索引（逻辑数组）
    validIdx = ismember(RAW_rxTime, NAV_rxTime);
    
    % 筛选出有效的RAW数据
    RAW_sorted = RAW(validIdx);
    
    % 为每个保留的数据计算txTime并添加为新字段
    for i = 1:length(RAW_sorted)
        RAW_sorted(i).txTime = RAW_sorted(i).rxTime - RAW_sorted(i).prRaw / c;
    end
    
end