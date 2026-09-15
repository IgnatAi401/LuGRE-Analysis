% clc; clear;
% load total.mat

n = numel(NAV_all);
x = 1:n;
y = [NAV_all.GDOP];

% 提取位置与 type
X = [NAV_all.posX];
Y = [NAV_all.posY];
Z = [NAV_all.posZ];
type = [NAV_all.type];

% 轨道高度（用于误差判断）
r = sqrt(X.^2 + Y.^2 + Z.^2);
idx_wrong = r >= 5e8;

% 按 type 分段（排除误差点）
idx_C = (type == 0) & ~idx_wrong;
idx_T = (type == 1) & ~idx_wrong;
idx_S = (type == 2) & ~idx_wrong;

figure(1);
hold on;
grid on;

scatter(x(idx_C), y(idx_C), 10, [0 0.4470 0.7410], 'filled');      % C
scatter(x(idx_T), y(idx_T), 10, [0.8500 0.3250 0.0980], 'filled');% T
scatter(x(idx_S), y(idx_S), 10, [0.4660 0.6740 0.1880], 'filled');% S
scatter(x(idx_wrong), y(idx_wrong), 15, [0 0 0], 'filled');      % Error

% xlabel('Epoch');
% ylabel('GDOP');
% title('GDOP vs Epoch (Segmented by Type, Error by Radius)');
ylim([0 100000]);

legend({'C', ...
        'T', ...
        'S', ...
        'Error'}, ...
        'Location', 'best');

y2 = [NAV_all.nSat];
x2 = 1:n;

figure(2);
hold on;
grid on;

scatter(x2(idx_C), y2(idx_C), 15, [0 0.4470 0.7410], 'filled');      % C
scatter(x2(idx_T), y2(idx_T), 15, [0.8500 0.3250 0.0980], 'filled');% T
scatter(x2(idx_S), y2(idx_S), 15, [0.4660 0.6740 0.1880], 'filled');% S
scatter(x2(idx_wrong), y2(idx_wrong), 20, [0 0 0], 'filled');      % Error

% xlabel('Epoch');
% ylabel('Number of Satellites');
% title('nSat vs Epoch (Segmented by Type, Error by Radius)');

legend({'C', ...
        'T', ...
        'S', ...
        'Error'}, ...
        'Location', 'best');
