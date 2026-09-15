function [r_ecef, v_ecef] = j2000_to_ecef(r_j2000, v_j2000, gpst_seconds)
    % 地球自转角速度 (rad/s) - WGS84 标准值
    omega_e = [0; 0; 7.2921151467e-5];
    
    % 获取旋转矩阵 (复用之前的逻辑)
    theta = calculate_gmst(gpst_seconds);
    R = [ cos(theta), sin(theta), 0;
         -sin(theta), cos(theta), 0;
          0,          0,          1];

    % 1. 位置转换
    r_ecef = R * r_j2000(:);
    
    % 2. 速度转换 (包含科里奥利项)
    % 公式: v_ecef = R * (v_j2000 - omega_e x r_j2000)
    v_ecef = R * (v_j2000(:) - cross(omega_e, r_j2000(:)));
end