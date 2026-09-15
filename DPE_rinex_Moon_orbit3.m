function [pos, dt, idx, dpos] = DPE_rinex_Moon_orbit3(satPos, Grid_ECEF, Origin_ECEF, Pred_Pos, k)
%==========================================================================
% Function: DPE for MOON Receiver (with Ionospheric Correction & Position Prior)
% Inputs:
%       satPos      - [n x 5] double
%       Grid_ECEF   - Candidate positions [4 x N_grid]
%       Origin_ECEF - Current estimated center [3 x 1]
%       Pred_Pos    - Predicted/Prior position [3 x 1] (NEW)
% Outputs:
%       pos, dt, idx, dpos
%==========================================================================
    % == Constants ===========================================================
    vel_c   = 299792458;          
    Omega_e = 7.2921151467e-5;    
    pi      = 3.141592653589793;
    
    % ---- Penalty Parameter (Hyperparameter) ----
    % lambda 代表预测位置的"硬度"。
    % 如果预测非常准（如来自惯导），lambda 调大；如果预测只是大概范围，调小。
    % 考虑到原代价函数中 10^(CN0/10) 量级通常在 10^3~10^5 之间，
    % lambda 设为平均 CN0 权重的 0.1~1 倍比较合理。
    lambda = 0.00001; 

    ion_Iz  = 5e-9;               
    n_sat   = size(satPos, 1);
    n_grid  = size(Grid_ECEF, 2);
    val_all = zeros(1, n_grid);

    % == 1. 计算卫星观测代价 (Original DPE Logic) =============================
    for i = 1 : n_sat
        curr_sat_pos = satPos(i, 1:3)';   
        curr_obs_pr  = satPos(i, 7);      
        cn0          = satPos(i, 8);      
        
        % Sagnac correction
        traveltime = curr_obs_pr / vel_c;
        angle = traveltime * Omega_e;
        R = [ cos(angle)  sin(angle) 0;
             -sin(angle)  cos(angle) 0;
              0           0          1 ];
        Rot_satpos = R * curr_sat_pos;
        
        % Elevation computation
        los = Rot_satpos - Origin_ECEF;
        r   = norm(los);
        el  = asin(los(3) / r) * 180 / pi;   
        
        % Ionospheric correction
        ion_F = 1 + 16 * (0.53 - el/180)^3;
        ion   = ion_Iz * ion_F * vel_c;      
        
        % Geometric distance & Model PR
        ri_candidate = vecnorm(Rot_satpos - Grid_ECEF(1:3, :));
        model_pr = ri_candidate + Grid_ECEF(4, :) + ion;
        
        % Cost accumulation (CN0 weighted L1 norm)
        % val_all = val_all + 10^(cn0/10) * abs(curr_obs_pr - model_pr);
        val_all = val_all + abs(curr_obs_pr - model_pr);
    end

    % == 2. 加入预测位置惩罚项 (New Penalty Logic) ============================
    % 计算网格点到预测位置的欧氏距离
    dist_to_pred = vecnorm(Grid_ECEF(1:3, :) - Pred_Pos);

    % 将惩罚项累加到代价函数中
    % 惩罚项 = lambda * 距离
    if k > 1
        val_all = val_all + lambda * dist_to_pred;
    end

    % == Minimum Search ======================================================
    [~, idx] = min(val_all);
    res  = Grid_ECEF(:, idx);
    pos  = res(1:3);
    dt   = res(4);
    dpos = pos - Origin_ECEF;
end