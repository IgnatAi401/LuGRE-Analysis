function [RAW_pos] = getPosRAW(RAW, gpsData, galileoData)
% GETPOSRAW 根据GPS总秒数匹配星历，计算卫星ECEF坐标和钟差(米)


    % --- 物理常数 ---
    C_LIGHT = 299792458;            % 光速 (m/s)
    OMEGA_E_DOT = 7.2921151467e-5;  % 地球自转角速度 (rad/s)
    
    % GPS 地心引力常数 (IS-GPS-200)
    MU_GPS = 3.986005e14; 
    
    % Galileo 地心引力常数 (Galileo OS SIS ICD)
    MU_GAL = 3.986004418e14; 

    % 初始化输出
    RAW_pos = RAW;
    numRaw = length(RAW);
    
    for i = 1:numRaw
        % 1. 获取当前观测数据
        t_tx_total = RAW(i).txTime; % 这里的 txTime 是 GPS总秒数
        sigId = RAW(i).signalId;
        prn = RAW(i).svId;
        
        targetSubset = [];
        mu_val = MU_GPS; % 默认常数
        
        % 2. 筛选对应星座和PRN的星历
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
        
        % 如果没有找到该星的任何星历
        if isempty(targetSubset)
            RAW_pos(i).ECEF = [NaN; NaN; NaN];
            RAW_pos(i).clockBias = NaN;
            disp('星历匹配出错')
            continue;
        end
        
        % 3. 核心匹配逻辑：基于 GPSTotalSeconds 找最近的星历
        % 假设星历结构体里的 GPSTotalSeconds 代表该星历参考时间(Toe)的绝对秒数
        [~, minIdx] = min(abs([targetSubset.GPSTotalSeconds] - t_tx_total));
        bestEph = targetSubset(minIdx);
        
        % 4. 计算 PVT (位置和钟差)
        % 传入总秒数进行计算
        [pos, vel, dts] = calculateSatPVT(bestEph, t_tx_total, mu_val, OMEGA_E_DOT, C_LIGHT);
        
        % 5. 保存结果
        RAW_pos(i).ECEF = pos;
        RAW_pos(i).Vel = vel;
        % 将钟差从秒转换为米
        RAW_pos(i).clockBias = dts; 
    end
end

