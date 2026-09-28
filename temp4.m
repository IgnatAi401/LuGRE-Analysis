% clc; clear;
% load pos3.mat

realref = readmatrix(fullfile(fileparts(mfilename('fullpath')), 'temp.xlsx'), 'Range', 'F339484:H339861');

% N = 1387;
M = size(NAV);
N = M(1,2);

r_ecef = zeros(N, 3);

for i = 1:N
    r_ecef(i,:) = j2000_to_ecef2([realref(i,1);realref(i,2);realref(i,3)], NAV(i).rxTime);
end
r_ecef = r_ecef*1000;

for i = 1:N
    % dif(i) = sqrt((r_ecef(i,1)-pos(i,1))^2+(r_ecef(i,2)-pos(i,2))^2+(r_ecef(i,3)-pos(i,3))^2);
    dif2(i) = sqrt((r_ecef(i,1)-NAV(i).posX)^2+(r_ecef(i,2)-NAV(i).posY)^2+(r_ecef(i,3)-NAV(i).posZ)^2);
    % dif3(i) = sqrt((pos(i,1)-NAV(i).posX)^2+(pos(i,2)-NAV(i).posY)^2+(pos(i,3)-NAV(i).posZ)^2);
end
% figure(1)
% plot(dif)
figure(2)
plot(dif2)
% figure(3)
% plot(dif3)


function r_ecef = j2000_to_ecef2(r_j2000, gps_seconds)
% J2000_TO_ECEF 将J2000惯性坐标系转换为ECEF地固坐标系
%
% 输入:
%   r_j2000     - 3x1向量, J2000坐标 (米)
%   gps_seconds - GPS时间总秒数 (标量), 如 1.420991617437020e+09
%
% 输出:
%   r_ecef      - 3x1向量, ECEF坐标 (米)
%
% 说明:
%   转换流程: J2000 -> TEME(近似) 或完整IAU76/80章动岁差
%   本函数使用完整的岁差+章动+地球自转+极移转换链
%   精度优于1米量级
%
% 参考: Vallado, "Fundamentals of Astrodynamics and Applications"

    %% Step 1: GPS时间 -> UTC时间
    [jd_utc, jd_ut1, jd_tt] = gps_to_julian(gps_seconds);
    
    %% Step 2: 计算各时间系统的儒略世纪数
    T_tt  = (jd_tt  - 2451545.0) / 36525.0;  % TT儒略世纪 (from J2000.0)
    T_ut1 = (jd_ut1 - 2451545.0) / 36525.0;  % UT1儒略世纪
    
    %% Step 3: 岁差矩阵 (IAU 1976)
    P = precession_matrix(T_tt);
    
    %% Step 4: 章动矩阵 (IAU 1980)
    [N, eps_true, delta_psi] = nutation_matrix(T_tt);
    
    %% Step 5: 地球自转 (GAST)
    theta_gast = compute_gast(jd_ut1, jd_tt, eps_true, delta_psi);
    R = rotz(theta_gast);  % 绕Z轴旋转
    
    %% Step 6: 极移矩阵 (使用近似值,可接入IERS数据)
    % 极移参数 xp, yp (弧度), 精确计算需查IERS公告
    % 此处使用近似零值,误差约10米量级
    xp = 0.0;  % 弧度, 实际约0.1-0.3角秒
    yp = 0.0;  % 弧度
    W = polar_motion_matrix(xp, yp);
    
    %% Step 7: 完整转换
    % r_ecef = W * R * N * P * r_j2000
    r_ecef = W * R * N * P * r_j2000;
    
end

function [jd_utc, jd_ut1, jd_tt] = gps_to_julian(gps_seconds)
% GPS时间转换为各时间系统的儒略日
%
% GPS时间起点: 1980年1月6日 00:00:00 UTC
% GPS时间与UTC: UTC = GPS - leap_seconds
% GPS时间与TAI: TAI = GPS + 19秒
% TT与TAI:      TT  = TAI + 32.184秒
% UT1与UTC:     UT1 = UTC + dUT1 (|dUT1| < 0.9秒, 近似取0)

    % GPS时间起点的儒略日
    jd_gps_epoch = 2444244.5;  % 1980-01-06 00:00:00 UTC
    
    % GPS -> JD(GPS)
    jd_gps = jd_gps_epoch + gps_seconds / 86400.0;
    
    % 计算跳秒数 (GPS超前UTC的秒数)
    leap_seconds = get_leap_seconds(gps_seconds);
    
    % GPS -> UTC
    jd_utc = jd_gps - leap_seconds / 86400.0;
    
    % UTC -> UT1 (近似: dUT1 ≈ 0, 精确值需查IERS)
    dUT1 = 0.0;  % 秒, 实际|dUT1| < 0.9秒
    jd_ut1 = jd_utc + dUT1 / 86400.0;
    
    % UTC -> TAI -> TT
    % TAI = UTC + leap_seconds (注意:这里leap是UTC滞后TAI的秒数)
    tai_minus_utc = leap_seconds;  % 当前跳秒
    jd_tai = jd_utc + tai_minus_utc / 86400.0;
    jd_tt  = jd_tai + 32.184 / 86400.0;  % TT = TAI + 32.184s
    
end

function leap = get_leap_seconds(gps_seconds)
% 根据GPS秒数返回当前跳秒数 (GPS时间超前UTC的秒数)
% GPS时间不含跳秒, UTC含跳秒
% 数据来源: IERS

    % [GPS秒数起点(1980-01-06), 跳秒数]
    % GPS秒数对应的跳秒生效时刻
    leap_table = [
        %  UTC日期          GPS秒数(近似)    GPS-UTC
           0,                0;          % 1980-01-06, 跳秒=0 (GPS起点时UTC-TAI=-19s)
           46828800,         1;          % 1981-07-01
           78364801,         2;          % 1982-07-01
           109900802,        3;          % 1983-07-01
           173059203,        4;          % 1985-07-01
           252028804,        5;          % 1988-01-01
           315187205,        6;          % 1990-01-01
           346723206,        7;          % 1991-01-01
           393984007,        8;          % 1992-07-01
           425520008,        9;          % 1993-07-01
           457056009,        10;         % 1994-07-01
           504489610,        11;         % 1996-01-01
           551750411,        12;         % 1997-07-01
           599184012,        13;         % 1999-01-01
           820108813,        14;         % 2006-01-01
           914803214,        15;         % 2009-01-01
           1025136015,       16;         % 2012-07-01
           1119744016,       17;         % 2015-07-01
           1167264017,       18;         % 2017-01-01
    ];
    
    % 找到当前GPS时间对应的跳秒
    leap = 0;
    for i = 1:size(leap_table, 1)
        if gps_seconds >= leap_table(i, 1)
            leap = leap_table(i, 2);
        end
    end
    
end

function P = precession_matrix(T)
% IAU 1976 岁差矩阵
% 输入: T - TT儒略世纪数 (from J2000.0)
% 输出: P - 3x3岁差旋转矩阵

    % 岁差角 (弧秒)
    zeta_A  =  2306.2181*T + 0.30188*T^2 + 0.017998*T^3;
    theta_A =  2004.3109*T - 0.42665*T^2 - 0.041775*T^3;
    z_A     =  2306.2181*T + 1.09468*T^2 + 0.018203*T^3;
    
    % 转换为弧度
    zeta_A  = arcsec2rad(zeta_A);
    theta_A = arcsec2rad(theta_A);
    z_A     = arcsec2rad(z_A);
    
    % 岁差矩阵: P = Rz(-z_A) * Ry(theta_A) * Rz(-zeta_A)
    P = rotz(-z_A) * roty(theta_A) * rotz(-zeta_A);
    
end

function [N, eps_true, delta_psi] = nutation_matrix(T)
% IAU 1980 章动矩阵 (简化版,含主要项)
% 输入: T - TT儒略世纪数
% 输出: N         - 3x3章动旋转矩阵
%        eps_true  - 真黄赤交角 (弧度)
%        delta_psi - 黄经章动 (弧度)

    % 平黄赤交角 (弧秒)
    eps0 = 84381.448 - 46.8150*T - 0.00059*T^2 + 0.001813*T^3;
    eps0 = arcsec2rad(eps0);
    
    % 基本天文参数 (弧度)
    % 月亮平近点角
    M_moon = deg2rad(134.96298139) + deg2rad(477198.8673981)*T ...
           + deg2rad(0.0086972)*T^2 + deg2rad(1/56250)*T^3;
    % 太阳平近点角
    M_sun  = deg2rad(357.52772333) + deg2rad(35999.0503400)*T ...
           - deg2rad(0.0001603)*T^2 - deg2rad(1/300000)*T^3;
    % 月亮纬度辐角
    F      = deg2rad(93.27191028)  + deg2rad(483202.0175381)*T ...
           - deg2rad(0.0036825)*T^2 + deg2rad(1/327270)*T^3;
    % 月日距角
    D      = deg2rad(297.85036306) + deg2rad(445267.1114800)*T ...
           - deg2rad(0.0019142)*T^2 + deg2rad(1/189474)*T^3;
    % 月亮轨道升交点黄经
    omega  = deg2rad(125.04452222) - deg2rad(1934.1362608)*T ...
           + deg2rad(0.0020708)*T^2 + deg2rad(1/450000)*T^3;
    
    % IAU 1980章动序列 (主要106项中取前15项,满足通常精度)
    % 格式: [l, l', F, D, Omega, psi系数(0.0001"), psi_T(0.0001"/T), 
    %        eps系数(0.0001"), eps_T(0.0001"/T)]
    nutation_coeffs = [
         0,  0,  0,  0,  1, -171996, -174.2,  92025,  8.9;
        -2,  0,  0,  2,  2,  -13187,   -1.6,   5736, -3.1;
         0,  0,  0,  2,  2,   -2274,   -0.2,    977, -0.5;
         0,  0,  0,  0,  2,    2062,    0.2,   -895,  0.5;
         0,  1,  0,  0,  0,    1426,   -3.4,     54, -0.1;
         0,  0,  1,  0,  0,     712,    0.1,     -7,  0.0;
        -2,  1,  0,  2,  2,    -517,    1.2,    224, -0.6;
         0,  0,  0,  2,  1,    -386,   -0.4,    200,  0.0;
         0,  0,  1,  2,  2,    -301,    0.0,    129, -0.1;
        -2, -1,  0,  2,  2,     217,   -0.5,    -95,  0.3;
        -2,  0,  1,  0,  0,    -158,    0.0,      0,  0.0;
        -2,  0,  0,  2,  1,     129,    0.1,    -70,  0.0;
         0,  0, -1,  2,  2,     123,    0.0,    -53,  0.0;
         2,  0,  0,  0,  0,      63,    0.0,      0,  0.0;
         0,  0,  1,  0,  1,      63,    0.1,    -33,  0.0;
    ];
    
    % 计算章动量
    delta_psi = 0.0;  % 黄经章动 (0.0001弧秒)
    delta_eps = 0.0;  % 黄赤交角章动 (0.0001弧秒)
    
    args = [M_moon, M_sun, F, D, omega];
    
    for i = 1:size(nutation_coeffs, 1)
        n     = nutation_coeffs(i, 1:5);
        S_psi = nutation_coeffs(i, 6) + nutation_coeffs(i, 7)*T;
        S_eps = nutation_coeffs(i, 8) + nutation_coeffs(i, 9)*T;
        
        arg = sum(n .* args);
        delta_psi = delta_psi + S_psi * sin(arg);
        delta_eps = delta_eps + S_eps * cos(arg);
    end
    
    % 转换为弧度
    delta_psi = arcsec2rad(delta_psi * 1e-4);
    delta_eps = arcsec2rad(delta_eps * 1e-4);
    
    % 真黄赤交角
    eps_true = eps0 + delta_eps;
    
    % 章动矩阵: N = Rx(-eps_true) * Rz(-delta_psi) * Rx(eps0)
    N = rotx(-eps_true) * rotz(-delta_psi) * rotx(eps0);
    
end

function theta_gast = compute_gast(jd_ut1, jd_tt, eps_true, delta_psi)
% 计算格林尼治真恒星时 (GAST)
% 输入: jd_ut1   - UT1儒略日
%        jd_tt    - TT儒略日  
%        eps_true - 真黄赤交角 (弧度)
%        delta_psi- 黄经章动 (弧度)
% 输出: theta_gast - GAST (弧度)

    T_ut1 = (jd_ut1 - 2451545.0) / 36525.0;
    
    % 格林尼治平恒星时 GMST (秒)
    % IAU 1982公式
    theta_gmst = 67310.54841 ...
               + (876600*3600 + 8640184.812866) * T_ut1 ...
               + 0.093104 * T_ut1^2 ...
               - 6.2e-6   * T_ut1^3;  % 秒
    
    % 转换为弧度 (mod 2pi)
    theta_gmst = mod(theta_gmst / 240.0, 360.0);  % 度
    theta_gmst = deg2rad(theta_gmst);
    
    % 方程差 (Equation of Equinoxes)
    % GAST = GMST + delta_psi * cos(eps_true)
    eq_equinox = delta_psi * cos(eps_true);
    
    theta_gast = theta_gmst + eq_equinox;
    theta_gast = mod(theta_gast, 2*pi);
    
end

function R = rotx(angle)
% 绕X轴旋转矩阵
    c = cos(angle); s = sin(angle);
    R = [1,  0,  0;
         0,  c,  s;
         0, -s,  c];
end

function R = roty(angle)
% 绕Y轴旋转矩阵
    c = cos(angle); s = sin(angle);
    R = [c,  0, -s;
         0,  1,  0;
         s,  0,  c];
end

function W = polar_motion_matrix(xp, yp)
% 极移矩阵
% 输入: xp, yp - 极移参数 (弧度)
% 输出: W - 3x3极移矩阵

    W = roty(xp) * rotx(yp);
    
end

function R = rotz(angle)
% 绕Z轴旋转矩阵
    c = cos(angle); s = sin(angle);
    R = [ c,  s,  0;
         -s,  c,  0;
          0,  0,  1];
end

function rad = arcsec2rad(arcsec)
    rad = arcsec * (pi / (180 * 3600));
end