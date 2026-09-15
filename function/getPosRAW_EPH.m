function [RAW_pos] = getPosRAW_EPH(RAW, gpsData, galileoData)
% GETPOSRAW_EPH
% 使用 RAW 中的 Toe，与星历中的 Toe 进行匹配，计算卫星 ECEF 坐标和钟差

    % --- 物理常数 ---
    C_LIGHT = 299792458;            % 光速 (m/s)
    OMEGA_E_DOT = 7.2921151467e-5;  % 地球自转角速度 (rad/s)
    
    % 地心引力常数
    MU_GPS = 3.986005e14;           % GPS
    MU_GAL = 3.986004418e14;        % Galileo

    % 初始化输出
    RAW_pos = RAW;
    numRaw = length(RAW);
    
    for i = 1:numRaw
        % 1. 当前观测
        toe_raw = RAW(i).Toe;     % ★ 使用 RAW 里的 Toe（周内秒）
        sigId   = RAW(i).signalId;
        prn     = RAW(i).svId;
        
        targetSubset = [];
        mu_val = MU_GPS; % 默认
        
        % 2. 根据星座筛选星历
        if (sigId == 0 || sigId == 1)
            % GPS
            if isfield(gpsData, 'PRN')
                targetSubset = gpsData([gpsData.PRN] == prn);
            end
            mu_val = MU_GPS;
            
        elseif (sigId == 2 || sigId == 3 || sigId == 4)
            % Galileo
            if isfield(galileoData, 'PRN')
                targetSubset = galileoData([galileoData.PRN] == prn);
            end
            mu_val = MU_GAL;
            
        else
            % 未知星座
            RAW_pos(i).ECEF = [NaN; NaN; NaN];
            RAW_pos(i).clockBias = NaN;
            continue;
        end
        
        % 3. 若该卫星没有任何星历
        if isempty(targetSubset)
            RAW_pos(i).ECEF = [NaN; NaN; NaN];
            RAW_pos(i).clockBias = NaN;
            disp('星历匹配出错：未找到对应 PRN')
            continue;
        end
        
        % 4. ★ 核心匹配逻辑：使用 Toe 匹配星历
        % 直接找 Toe 最接近的那条星历
        ephToe = [targetSubset.Toe];
        [~, minIdx] = min(abs(ephToe - toe_raw));
        bestEph = targetSubset(minIdx);
        
        % 5. 计算卫星位置与钟差
        % 注意：这里仍然用 txTime（GPS 总秒数）进行轨道传播
        t_tx_total = RAW(i).txTime;
        [pos, vel, dts] = calculateSatPVT( ...
            bestEph, t_tx_total, mu_val, OMEGA_E_DOT, C_LIGHT);
        
        % 6. 保存结果
        RAW_pos(i).ECEF = pos;
        RAW_pos(i).Vel = vel;
        RAW_pos(i).clockBias = dts;  % 秒或米取决于 calculateSatPVT 的实现
    end
end
