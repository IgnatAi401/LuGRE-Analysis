function satPos = getPos(rxTime, RAW)
% GETPOS 根据 rxTime 匹配数据，输出卫星坐标和修正后的伪距
%
% 输入:
%   rxTime - 查询时间 (double)
%   RAW    - 1xn struct，包含字段: rxTime, ECEF, prRaw, clockBias
%
% 输出:
%   satPos - Nx4 矩阵
%            前3列: [X, Y, Z] (ECEF坐标)
%            第4列: [prRaw + clockBias] (修正后的伪距)

    % 1. 提取所有时间并查找匹配索引
    allTimes = [RAW.rxTime];
    
    % 使用微小误差范围判断浮点时间是否相等
    tol = 1e-12;
    idx = abs(allTimes - rxTime) < tol;
    
    % 2. 提取并处理数据
    if ~any(idx)
        satPos = []; % 未找到匹配数据
    else
        % --- 处理 ECEF (前3列) ---
        % 将匹配的ECEF数据拼接后重塑为 Nx3 矩阵
        ecefTemp = [RAW(idx).ECEF];
        ecefMat = reshape(ecefTemp, 3, []).'; 

        velTemp = [RAW(idx).Vel];
        velMat = reshape(velTemp, 3, []).'; 
        
        % --- 处理 第4列 (prRaw + clockBias) ---
        % 提取伪距
        prVec = [RAW(idx).prRaw];
        % 提取钟差
        biasVec = [RAW(idx).clockBias];
        % 提取cn0
        cn0vec = [RAW(idx).cn0]';

        fdvec = [RAW(idx).fdRaw]';

        fdratevec = [RAW(idx).fdRateRaw]';

        cpvec = [RAW(idx).carrierPhase]';
        
        % 计算修正后的伪距 (确保转为列向量相加)
        correctedPr = prVec(:) + biasVec(:);
        
        % 3. 拼接最终结果 [X, Y, Z, (prRaw + clockBias)]
        satPos_Raw = [ecefMat, velMat, correctedPr, cn0vec, fdvec, fdratevec, cpvec];
        [~, idx] = unique(satPos_Raw(:,4), 'stable');
        satPos = satPos_Raw(idx, :);

    end

end