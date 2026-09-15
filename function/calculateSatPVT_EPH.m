function [pos, dts] = calculateSatPVT_EPH(eph, t_total_current, mu, omega_e_dot, c)
% CALCULATE SAT PVT
% 输入:
%   eph: 选中的星历结构体
%   t_total_current: 当前信号发射时刻 (GPS总秒数)
%   mu, omega_e_dot, c: 物理常数
% 输出:
%   pos: ECEF坐标 (m)
%   dts: 卫星钟差 (m)

    % --- 参数提取 ---
    sqrtA   = eph.SqrtA;
    e       = eph.e;
    %这里不再使用 eph.Toe 来计算时间差，而是使用 GPSTotalSeconds
    t_total_oe = eph.Toe; 
    
    M0      = eph.M0;
    delta_n = eph.DeltaN;
    
    i0      = eph.i0;
    omega0  = eph.Omega0;
    omega   = eph.omega;
    dot_i   = eph.IDOT;
    dot_omega = eph.OmegaDot; 
    
    cuc = eph.Cuc; cus = eph.Cus;
    crc = eph.Crc; crs = eph.Crs;
    cic = eph.Cic; cis = eph.Cis;
    
    af0 = eph.ClockBias;
    af1 = eph.ClockDrift;
    af2 = eph.ClockDriftRate;
    
    % --- 1. 时间计算 ---
    % 计算 t_k
    tk = t_total_current - t_total_oe;
   
    
    % 理论上如果匹配正确，tk应该很小（+/- 2小时内）。
    % 考虑到极端的匹配情况或数据异常，仍保留 +/- 302400 的归一化逻辑作为保险，
    % if tk > 302400
    %     tk = tk - 604800;
    % elseif tk < -302400
    %     tk = tk + 604800;
    % end
    
    % --- 2. 轨道平运动计算 ---
    A = sqrtA^2;
    n0 = sqrt(mu / A^3);
    n = n0 + delta_n;
    
    % 平近点角
    Mk = M0 + n * tk;
    
    % --- 3. 开普勒方程迭代 (高精度) ---
    Ek = Mk; 
    for k = 1:20
        Ek_old = Ek;
        Ek = Mk + e * sin(Ek);
        if abs(Ek - Ek_old) < 1e-14
            break;
        end
    end
    
    % --- 4. 钟差计算 ---
    % 相对论效应常数 F
    F = -2 * sqrt(mu) / c^2; 
    dtr = F * e * sqrtA * sin(Ek);
    
    % 总钟差 (米)
    dts = (af0 + af1 * tk + af2 * tk^2 + dtr)*299792458;
    
    % --- 5. 坐标计算 ---
    % 真近点角
    sin_v = (sqrt(1 - e^2) * sin(Ek)) / (1 - e * cos(Ek));
    cos_v = (cos(Ek) - e) / (1 - e * cos(Ek));
    vk = atan2(sin_v, cos_v);
    
    % 升交点角距
    phi_k = vk + omega;
    
    % 二阶谐波摄动
    sin_2phi = sin(2 * phi_k);
    cos_2phi = cos(2 * phi_k);
    
    du_k = cus * sin_2phi + cuc * cos_2phi;
    dr_k = crs * sin_2phi + crc * cos_2phi;
    di_k = cis * sin_2phi + cic * cos_2phi;
    
    % 校正后的参数
    uk = phi_k + du_k;
    rk = A * (1 - e * cos(Ek)) + dr_k;
    ik = i0 + dot_i * tk + di_k;
    
    % 轨道平面位置
    xk_prime = rk * cos(uk);
    yk_prime = rk * sin(uk);
    
    % 升交点经度 (需要用到 Toe 的周内秒部分进行地球自转校正)
    
    toe_week_sec = eph.Toe; 
    Omega_k = omega0 + (dot_omega - omega_e_dot) * tk - omega_e_dot * toe_week_sec;
    
    % ECEF 坐标转换
    x = xk_prime * cos(Omega_k) - yk_prime * cos(ik) * sin(Omega_k);
    y = xk_prime * sin(Omega_k) + yk_prime * cos(ik) * cos(Omega_k);
    z = yk_prime * sin(ik);
    
    pos = [x; y; z];
end