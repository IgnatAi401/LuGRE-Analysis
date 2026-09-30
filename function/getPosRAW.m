function [RAW_pos] = getPosRAW(RAW, gpsData, galileoData, maxAgeS)
% GETPOSRAW 根据GPS总秒数匹配星历，计算卫星ECEF坐标和钟差(米)
% maxAgeS：发射时刻与 Toe 的最大允许间隔(s)，默认 4 h；超出置 NaN，
% 防止误用其他日期的星历时静默得到 km 级错误的卫星位置。
    if nargin < 4, maxAgeS = 4*3600; end
    missingCount = 0; staleCount = 0;


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
        RAW_pos(i).ECEF = nan(3,1);
        RAW_pos(i).Vel = nan(3,1);
        RAW_pos(i).clockBias = NaN;
        t_tx_total = RAW(i).txTime; % 这里的 txTime 是 GPS总秒数
        sigId = RAW(i).signalId;
        prn = RAW(i).svId;
        
        targetSubset = [];
        
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
            missingCount = missingCount + 1;
            continue;
        end

        % 3. 基于轨道参考时刻 Toe 找最近的星历
        % GPSTotalSeconds 是 Toc；匹配实际轨道参考时刻 Toe。
        toc = [targetSubset.GPSTotalSeconds];
        toe = [targetSubset.Toe] + 604800*round((toc-[targetSubset.Toe])/604800);
        [age, minIdx] = min(abs(toe - t_tx_total));
        if age > maxAgeS
            staleCount = staleCount + 1;
            continue;
        end
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
    if missingCount > 0
        warning('LuGRE:EphemerisMissing','%d 条观测找不到对应卫星的星历，已置 NaN。',missingCount);
    end
    if staleCount > 0
        warning('LuGRE:EphemerisStale', ...
            '%d/%d 条观测的最近星历距发射时刻超过 %.1f h，已置 NaN；请核对 RINEX 日期。', ...
            staleCount, numRaw, maxAgeS/3600);
    end
end

