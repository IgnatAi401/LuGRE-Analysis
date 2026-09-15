function [pos, vel, dts] = calculateSatPVT(eph, t_total_current, mu, omega_e_dot, c)
% CALCULATE SAT PVT
% 输入:
%   eph: 选中的星历结构体
%   t_total_current: 当前信号发射时刻 (GPS总秒数)
%   mu, omega_e_dot, c: 物理常数
% 输出:
%   pos: ECEF坐标 (m)
%   vel: ECEF速度 (m/s)
%   dts: 卫星钟差 (m)

    % --- 参数提取 ---
    sqrtA   = eph.SqrtA;
    e       = eph.e;
    t_total_oe = eph.GPSTotalSeconds; 
    
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

    % --- 时间 ---
    tk = t_total_current - t_total_oe;

    % --- 轨道平运动 ---
    A = sqrtA^2;
    n0 = sqrt(mu / A^3);
    n = n0 + delta_n;

    Mk = M0 + n * tk;

    % --- 开普勒方程 ---
    Ek = Mk;
    for k = 1:20
        Ek_old = Ek;
        Ek = Mk + e * sin(Ek);
        if abs(Ek - Ek_old) < 1e-14
            break;
        end
    end

    % --- 钟差 ---
    F = -2 * sqrt(mu) / c^2; 
    dtr = F * e * sqrtA * sin(Ek);
    dts = (af0 + af1 * tk + af2 * tk^2 + dtr)*299792458;

    % ======================
    % 位置计算
    % ======================

    sin_v = (sqrt(1 - e^2) * sin(Ek)) / (1 - e * cos(Ek));
    cos_v = (cos(Ek) - e) / (1 - e * cos(Ek));
    vk = atan2(sin_v, cos_v);

    phi_k = vk + omega;

    sin_2phi = sin(2 * phi_k);
    cos_2phi = cos(2 * phi_k);

    du_k = cus * sin_2phi + cuc * cos_2phi;
    dr_k = crs * sin_2phi + crc * cos_2phi;
    di_k = cis * sin_2phi + cic * cos_2phi;

    uk = phi_k + du_k;
    rk = A * (1 - e * cos(Ek)) + dr_k;
    ik = i0 + dot_i * tk + di_k;

    xk_prime = rk * cos(uk);
    yk_prime = rk * sin(uk);

    toe_week_sec = eph.Toe; 
    Omega_k = omega0 + (dot_omega - omega_e_dot) * tk - omega_e_dot * toe_week_sec;

    x = xk_prime * cos(Omega_k) - yk_prime * cos(ik) * sin(Omega_k);
    y = xk_prime * sin(Omega_k) + yk_prime * cos(ik) * cos(Omega_k);
    z = yk_prime * sin(ik);

    pos = [x; y; z];

    % ======================
    % 速度计算
    % ======================

    % 偏近点角速度
    Ek_dot = n / (1 - e*cos(Ek));

    % 真近点角速度
    vk_dot = sin(Ek)*Ek_dot*(1+e*cos(vk))/(sin(vk)*(1-e*cos(Ek)));

    % 摄动导数
    du_dot = 2*(cus*cos_2phi - cuc*sin_2phi)*vk_dot;
    dr_dot = 2*(crs*cos_2phi - crc*sin_2phi)*vk_dot;
    di_dot = 2*(cis*cos_2phi - cic*sin_2phi)*vk_dot;

    % 参数导数
    u_dot = vk_dot + du_dot;
    r_dot = A*e*sin(Ek)*Ek_dot + dr_dot;
    i_dot = dot_i + di_dot;

    % 轨道平面速度
    xk_dot = r_dot*cos(uk) - rk*sin(uk)*u_dot;
    yk_dot = r_dot*sin(uk) + rk*cos(uk)*u_dot;

    Omega_dot_k = dot_omega - omega_e_dot;

    % ECEF速度
    vx = xk_dot*cos(Omega_k) ...
       - yk_dot*cos(ik)*sin(Omega_k) ...
       + yk_prime*sin(ik)*sin(Omega_k)*i_dot ...
       - (xk_prime*sin(Omega_k) + yk_prime*cos(ik)*cos(Omega_k))*Omega_dot_k;

    vy = xk_dot*sin(Omega_k) ...
       + yk_dot*cos(ik)*cos(Omega_k) ...
       - yk_prime*sin(ik)*cos(Omega_k)*i_dot ...
       + (xk_prime*cos(Omega_k) - yk_prime*cos(ik)*sin(Omega_k))*Omega_dot_k;

    vz = yk_dot*sin(ik) + yk_prime*cos(ik)*i_dot;

    vel = [vx; vy; vz];

end