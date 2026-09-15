function gmst_rad = calculate_gmst(gpst_total_seconds)
    % 输入: 自 1980-01-06 00:00:00 以来累积的 GPS 秒数
    % 输出: 格林尼治平恒星时 (弧度)

    % 1. 将 GPS 秒转换为 UTC/UT1 下的儒略日期 (JD)
    % GPS 历元 JD 为 2444244.5。
    % 减去 18 秒跳秒 (GPS 比 UTC 快 18 秒)，将原子时对齐到世界时。
    leap_seconds = 18; 
    jd_ut1 = 2444244.5 + (gpst_total_seconds - leap_seconds) / 86400;

    % 2. 计算自 J2000.0 起算的世纪数 T
    T = (jd_ut1 - 2451545.0) / 36525.0;

    % 3. 计算 GMST (秒为单位的算法)
    gmst_sec = 67310.54841 + (876600 * 3600 + 8640184.812866) * T ...
               + 0.093104 * T^2 - 6.2e-6 * T^3;
    
    % 4. 转换为 [0, 2*pi] 弧度
    gmst_deg = mod(gmst_sec / 240, 360);
    gmst_rad = deg2rad(gmst_deg);
end