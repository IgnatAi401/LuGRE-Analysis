%clc; clear;
%load total.mat

n = numel(NAV_all);
X2   = zeros(n,1);
Y2   = zeros(n,1);
Z2   = zeros(n,1);
type = zeros(n,1);

for i = 1:n
    X2(i)   = NAV_all(i).posX;
    Y2(i)   = NAV_all(i).posY;
    Z2(i)   = NAV_all(i).posZ;
    type(i) = NAV_all(i).type;
end

figure;
hold on;
grid on;
axis equal;

% =====================
% 计算地心距（用于误差判断）
% =====================
r = sqrt(X2.^2 + Y2.^2 + Z2.^2);

% =====================
% 误差点（仍按原轨道高度逻辑）
% =====================
idx_wrong = r >= 5e8;

% =====================
% 按 type 划分阶段（不包含误差点）
% =====================
idx_C = (type == 0) & ~idx_wrong;
idx_T = (type == 1) & ~idx_wrong;
idx_S = (type == 2) & ~idx_wrong;

% =====================
% 分阶段绘制
% =====================
hC = scatter3(X2(idx_C), Y2(idx_C), Z2(idx_C), ...
              30, [0 0.4470 0.7410], 'filled');      % C

hT = scatter3(X2(idx_T), Y2(idx_T), Z2(idx_T), ...
              30, [0.8500 0.3250 0.0980], 'filled'); % T

hS = scatter3(X2(idx_S), Y2(idx_S), Z2(idx_S), ...
              30, [0.4660 0.6740 0.1880], 'filled'); % S

% =====================
% 误差点绘制（黑色）
% =====================
hW = scatter3(X2(idx_wrong), Y2(idx_wrong), Z2(idx_wrong), ...
              40, [0 0 0], 'filled');

% =====================
% 地心
% =====================
hE = scatter3(0, 0, 0, 80, 'r', 'filled');

% =====================
% 坐标与标题
% =====================
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
title('Earth–Moon Transfer Trajectory (Segmented by Type, Error by Radius)');

legend([hC hT hS hE hW], ...
       {'C phase (type=0)', ...
        'T phase (type=1)', ...
        'S phase (type=2)', ...
        'Earth center', ...
        'Error (r \ge 5e8 m)'}, ...
       'Location', 'best');

% legend([hC hT hS hE], ...
%        {'C phase (type=0)', ...
%         'T phase (type=1)', ...
%         'S phase (type=2)', ...
%         'Earth center'}, ...
%        'Location', 'best');

view(3);

% =====================
% 输出误差点历元号
% =====================
error_epochs = find(idx_wrong);

fprintf('发现 %d 个误差点（按轨道高度判定），历元号为：\n', numel(error_epochs));
disp(error_epochs.');
