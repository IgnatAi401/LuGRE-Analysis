% n=60;
% n1 = 1387 -n;
% n2 = n1 +1;
% t_gps_array = zeros(1,n);
% for i = 1:n
%     t_gps_array(i) = NAV(i+n1).rxTime;
% end
% tspan = NAV(1387).rxTime - NAV(n2).rxTime;
% moonEph = getMoonEph(t_gps_array)*1000;
% 
% [Xr,Vr] = ecef_to_j2000([NAV(n2).posX;NAV(n2).posY;NAV(n2).posZ],[NAV(n2).velX;NAV(n2).velY;NAV(n2).velZ],NAV(n2).rxTime);
% 
% X0 = [Xr; Vr; moonEph(1,2); moonEph(1,3); moonEph(1,4); moonEph(1,5); moonEph(1,6); moonEph(1,7)];
% [t, X] = ode45(@(t,X) Dynamic_R3BJ2(t,X,[]), [0 tspan], X0);
% % [t, X] = ode45(@Dynamic_R3BJ2, [0 1], []);
% nx = size(X,1);
% X2 = zeros(n,1);
% Y2 = zeros(n,1);
% Z2 = zeros(n,1);
% X3 = zeros(nx,1);
% Y3 = zeros(nx,1);
% Z3 = zeros(nx,1);
% for i = 1:n
%     X2(i) = NAV(i+n1).posX;
%     Y2(i) = NAV(i+n1).posY;
%     Z2(i) = NAV(i+n1).posZ;
% end
% 
% for i = 1:nx
%     [Xr, ~] = j2000_to_ecef([X(i,1);X(i,2);X(i,3)],[0;0;0],t(i)+NAV(n2).rxTime);
%     X3(i) = Xr(1);
%     Y3(i) = Xr(2);
%     Z3(i) = Xr(3);
% end
% figure;
% hold on;
% grid on;
% axis equal;
% 
% scatter3(X2, Y2, Z2, ...
%          10, 'b', 'filled');
% 
% scatter3(X3, Y3, Z3, ...
%          40, 'r', 'filled');
% view(3);
% 
% % 假设数据格式：
% % traj1: [N1 x 3] -> (X, Y, Z)
% % traj2: [N2 x 3] -> (X, Y, Z)
% tra1 = [X2,Y2,Z2];
% tra2 = [X3,Y3,Z3];
% [rmse, error_dist] = evaluate_trajectory_error(tra1, tra2);
% 
% function [rmse, error_dist] = evaluate_trajectory_error(traj1, traj2)
%     % 获取两组轨迹的点数
%     N1 = size(traj1, 1);
%     N2 = size(traj2, 1);
% 
%     % 1. 建立归一化的弧长（或者时间轴）作为插值基准
%     % 假设起点到终点的总进度为 0 到 1
%     t1 = linspace(0, 1, N1)';
%     t2 = linspace(0, 1, N2)';
% 
%     % 2. 将 traj2 插值到 traj1 的采样点上
%     % 使用 spline (样条插值) 平滑度更高，也可选 'linear'
%     traj2_interp = zeros(N1, 3);
%     traj2_interp(:,1) = interp1(t2, traj2(:,1), t1, 'spline');
%     traj2_interp(:,2) = interp1(t2, traj2(:,2), t1, 'spline');
%     traj2_interp(:,3) = interp1(t2, traj2(:,3), t1, 'spline');
% 
%     % 3. 计算逐点欧氏距离 (Error Distance)
%     % formula: sqrt((x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2)
%     diff = traj1 - traj2_interp;
%     error_dist = sqrt(sum(diff.^2, 2));
% 
%     % 4. 计算指标
%     rmse = sqrt(mean(error_dist.^2));
%     mae = mean(error_dist);
%     max_err = max(error_dist);
% 
%     % 打印结果
%     fprintf('--- 轨迹误差评估 ---\n');
%     fprintf('RMSE: %.4f\n', rmse);
%     fprintf('MAE:  %.4f\n', mae);
%     fprintf('Max Error: %.4f\n', max_err);
% 
%     % 5. 可视化
%     figure;
%     plot3(traj1(:,1), traj1(:,2), traj1(:,3), 'b-', 'LineWidth', 1.5); hold on;
%     plot3(traj2(:,1), traj2(:,2), traj2(:,3), 'r--', 'LineWidth', 1);
%     grid on; legend('轨迹 1 (基准)', '轨迹 2 (对比)');
%     title('轨迹空间对比');
% end
for i = 1:1387
    dx3(i) = abs(pos2(i,1) - NAV(i).posX);
    dy3(i) = abs(pos2(i,2) - NAV(i).posY);
    dz3(i) = abs(pos2(i,3) - NAV(i).posZ);
    dd3(i) = sqrt(dx3(i)^2+dy3(i)^2+dz3(i)^2);
end
figure(5)
plot(1:n, dd3)
xlabel('Epoch')
ylabel('3D Positioning Error / m')
title('LS Residuals Relative to Receiver Algorithm Solution')