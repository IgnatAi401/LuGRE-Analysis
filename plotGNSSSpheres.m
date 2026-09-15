function plotGNSSSpheresAndPlane(dataStruct, refPos, idx_start, idx_end)
% plotGNSSSpheresAndPlane 
% 功能1：绘制3D ECEF球壳
% 功能2：绘制参考高度平面上的投影圆（位置线）
%
% 输入:
%   dataStruct - 1*n struct, 包含 .ECEF (3x1) 和 .prRaw (1x1)
%   refPos     - 参考坐标 3x1 ECEF
%   idx_start, idx_end - 索引范围

    % 提取选定的卫星数量
    indices = idx_start:idx_end;
    numSat = length(indices);
    colors = lines(numSat);

    % 创建图形窗口，使用分栏布局
    figure('Name', 'GNSS Visualization', 'Color', 'w', 'Position', [100, 100, 1200, 600]);
    t = tiledlayout(1, 2, 'TileSpacing', 'compact');

    %% --- 左图：3D ECEF 视图 ---
    ax1 = nexttile;
    hold(ax1, 'on'); grid(ax1, 'on'); axis(ax1, 'equal');
    xlabel(ax1, 'ECEF X (m)'); ylabel(ax1, 'ECEF Y (m)'); zlabel(ax1, 'ECEF Z (m)');
    title(ax1, '3D 视图: 卫星定位球壳');
    view(ax1, 3);

    % 生成单位球数据
    [unitX, unitY, unitZ] = sphere(30);

    %% --- 右图：2D 平面投影视图 (Local Tangent Plane) ---
    ax2 = nexttile;
    hold(ax2, 'on'); grid(ax2, 'on'); axis(ax2, 'equal');
    xlabel(ax2, 'East (m)'); ylabel(ax2, 'North (m)');
    title(ax2, '2D 平面图: 伪距球壳与参考平面的交线');
    % 标出参考点 (原点)
    plot(ax2, 0, 0, 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r', 'DisplayName', '参考点');

    %% --- 循环处理数据 ---
    legendEntries = {}; % 用于存储图例
    
    count = 1;
    for k = indices
        if k > length(dataStruct), continue; end
        
        % 1. 获取数据
        satECEF = dataStruct(k).ECEF;
        radius = dataStruct(k).prRaw;
        col = colors(count, :);
        
        % ----------------------
        % 绘制 3D 球壳 (左图)
        % ----------------------
        X = unitX * radius + satECEF(1);
        Y = unitY * radius + satECEF(2);
        Z = unitZ * radius + satECEF(3);
        
        surf(ax1, X, Y, Z, 'FaceColor', col, 'FaceAlpha', 0.1, 'EdgeColor', 'none');
        plot3(ax1, satECEF(1), satECEF(2), satECEF(3), '.', 'Color', col, 'MarkerSize', 12);
        
        % ----------------------
        % 计算并绘制 2D 投影 (右图)
        % ----------------------
        % 步骤 A: 将卫星 ECEF 转换到以 refPos 为原点的 ENU 坐标系
        [E_sat, N_sat, U_sat] = ecef2enu_custom(satECEF, refPos);
        
        % 步骤 B: 计算截面圆的半径
        % 勾股定理: r_2d^2 + U_sat^2 = radius^2
        % r_2d = sqrt(radius^2 - U_sat^2)
        % U_sat 是卫星相对于参考平面的垂直高度
        
        if radius > abs(U_sat)
            r_2d = sqrt(radius^2 - U_sat^2);
            
            % 步骤 C: 在 2D 图上画圆
            % 圆心是卫星在平面的投影 (E_sat, N_sat)，半径是 r_2d
            theta = linspace(0, 2*pi, 200);
            circ_x = E_sat + r_2d * cos(theta);
            circ_y = N_sat + r_2d * sin(theta);
            
            % 绘制圆弧
            plot(ax2, circ_x, circ_y, 'Color', col, 'LineWidth', 1.5, ...
                'DisplayName', sprintf('Sat %d', k));
            
            % 同时也画一条线连向卫星方向(可选，指示卫星方位)
            % quiver(ax2, 0, 0, E_sat, N_sat, 'Color', [0.8 0.8 0.8], 'MaxHeadSize', 0);
        else
            warning('卫星 %d 的伪距小于其高度，无法与平面相交 (数据可能异常)', k);
        end
        
        count = count + 1;
    end

    %% --- 后期修饰 ---
    
    % 3D图修饰：画参考点
    plot3(ax1, refPos(1), refPos(2), refPos(3), 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
    
    % 2D图修饰：设置合适的视野
    % 因为卫星离得很远，圆非常大。为了看清交点，我们需要限制视野在参考点附近
    xlim(ax2, [-300, 300]); % 只看参考点周围 +/- 300米范围
    ylim(ax2, [-300, 300]);
    text(ax2, 0, 0, '  Ref', 'Color', 'r', 'FontWeight', 'bold');
    legend(ax2, 'Location', 'bestoutside');
    
    fprintf('绘图完成。\n注意：右侧2D图中，已自动缩放至参考点附近 +/-300米，\n以便观察各圆弧是否交汇于原点。\n');
end

%% --- 辅助函数：ECEF 转 ENU (不依赖工具箱) ---
function [E, N, U] = ecef2enu_custom(satECEF, refECEF)
    % 输入均为 3x1 列向量
    % 1. 计算参考点的经纬度 (用于构建旋转矩阵)
    dx = satECEF - refECEF;
    
    X = refECEF(1); Y = refECEF(2); Z = refECEF(3);
    
    % 简单的球体近似计算经纬度 (对于旋转矩阵精度足够)
    lam = atan2(Y, X);          % 经度
    phi = atan2(Z, sqrt(X^2 + Y^2)); % 纬度
    
    % 2. 旋转矩阵 (ECEF -> ENU)
    % R = [ -sin(lam)           cos(lam)          0     ;
    %       -sin(phi)cos(lam)  -sin(phi)sin(lam)  cos(phi);
    %        cos(phi)cos(lam)   cos(phi)sin(lam)  sin(phi)];
    
    sinL = sin(lam); cosL = cos(lam);
    sinP = sin(phi); cosP = cos(phi);
    
    R = [-sinL,        cosL,       0;
         -sinP*cosL,  -sinP*sinL,  cosP;
          cosP*cosL,   cosP*sinL,  sinP];
      
    % 3. 旋转差分向量
    enu = R * dx;
    E = enu(1);
    N = enu(2);
    U = enu(3);
end