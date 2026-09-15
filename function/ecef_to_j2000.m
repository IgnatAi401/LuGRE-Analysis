function [r_j2000, v_j2000] = ecef_to_j2000(r_ecef, v_ecef, gpst_seconds)
    omega_e = [0; 0; 7.2921151467e-5];
    
    theta = calculate_gmst(gpst_seconds);
    R = [ cos(theta), sin(theta), 0;
         -sin(theta), cos(theta), 0;
          0,          0,          1];

    % 1. 位置转换 (R 的转置)
    r_j2000 = R' * r_ecef(:);
    
    % 2. 速度转换
    % 公式: v_j2000 = R' * v_ecef + omega_e x r_j2000
    v_j2000 = R' * v_ecef(:) + cross(omega_e, r_j2000);
end