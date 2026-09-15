clc;clear;

addpath function
addpath reader
addpath settings

[ACQfilePath, RAWfilePath, NAVfilePath] = FindTelemetryTxtFiles('/Users/1gnat4y/Downloads/LuGRE', 'OP1_0');
disp(ACQfilePath)
disp(RAWfilePath)
disp(NAVfilePath)
[ACQ] = TxtParserAcq(ACQfilePath);
[RAW] = TxtParserRaw(RAWfilePath);
[NAV] = TxtParserNav(NAVfilePath);


% load matlab3.mat

[RAW_sorted] = sortRAW(NAV, RAW);



% csv_filepath = '/Users/1gnat4y/Downloads/LuGRE/L0/TLM/TLM_EPH_20250115_152133_01H_C_OP1_0.csv';
% RAW_sorted = match_toe(RAW_sorted, csv_filepath);



filename = "BRDC00IGS_R_20250150000_01D_MN.rnx";
[gpsData, galileoData] = read_gnss_rinex(filename);
[RAW_Pos] = getPosRAW(RAW_sorted, gpsData, galileoData);
% [RAW_Pos] = getPosRAW_EPH(RAW_sorted, gpsData, galileoData);

n = numel(NAV);

%n = 500;


t_gps_array = zeros(1,n);
for i = 1:n
    t_gps_array(i) = NAV(i).rxTime;
end

moonEph = getMoonEph(t_gps_array)*1000;







% n = 500;

dx = nan(n,1);
dy = nan(n,1);
dz = nan(n,1);
dx2 = nan(n,1);
dy2 = nan(n,1);
dz2 = nan(n,1);
dx3 = nan(n,1);
dy3 = nan(n,1);
dz3 = nan(n,1);
dt = nan(n,1);
dt2 = nan(n,1);
dt3 = nan(n,1);
dd = nan(n,1);
dd2 = nan(n,1);
dd3 = nan(n,1);
pos =nan(n,3);
pos2 =nan(n,3);
pos3 =nan(n,3);

% i =41;
%k =40;
% Origin_ECEF = [NAV(50).posX; NAV(50).posY; NAV(50).posZ];
Origin_ECEF = [NAV(1).posX; NAV(1).posY; NAV(1).posZ];
disp(norm(Origin_ECEF))

dt0 = -140000;
% pos0 = [-10000;-350000;-30000];
% Origin_ECEF = Origin_ECEF + pos0;

fprintf('Processing ...');
hwait = waitbar(0,'1','name','Calculation');
[Xr,Vr] = ecef_to_j2000([NAV(181).posX;NAV(181).posY;NAV(181).posZ],[NAV(181).velX;NAV(181).velY;NAV(181).velZ],NAV(181).rxTime);

X0 = [Xr; Vr; moonEph(4,2); moonEph(4,3); moonEph(4,4); moonEph(4,5); moonEph(4,6); moonEph(4,7)];

for i = 1:(n-1)
    satPos = getPos(NAV(i).rxTime, RAW_Pos);
    satVel = getVel(NAV(i).rxTime, RAW_Pos);



    Grid_ECEF = Grid_Generation_ECEF2(Origin_ECEF, i);
    Grid_ECEF (4,:) = Grid_ECEF (4,:) + dt0;

    count_num1 = floor((i-2)/60)*60 + 1;
    count_num2 = floor((i-1)/60)*60 + 1;
    if count_num1 < 0
        count_num1 = 1;
    end

    

    % if i > 201
    %     %% 4. 调用 ODE
    %     tspan = NAV(i).rxTime - NAV(201).rxTime;
    %     [t, X] = ode45(@(t,X) Dynamic_R3BJ2(t,X,[]), [0 tspan], X0);
    %     num_X = size(X,1);
    %     [relp, velp] = j2000_to_ecef([X(num_X,1);X(num_X,2);X(num_X,3)], [X(num_X,4);X(num_X,5);X(num_X,6)], NAV(i).rxTime);
    % else
    %     relp = [NAV(1).posX;NAV(1).posY;NAV(1).posZ];
    %     velp = [NAV(1).velX;NAV(1).velY;NAV(1).velZ];
    % end

    if i > 1
        %% 4. 调用 ODE
        tspan = NAV(i).rxTime - NAV(count_num1).rxTime;
        [t, X] = ode45(@(t,X) Dynamic_R3BJ2(t,X,[]), [0 tspan], X0);
        num_X = size(X,1);
        [relp, velp] = j2000_to_ecef([X(num_X,1);X(num_X,2);X(num_X,3)], [X(num_X,4);X(num_X,5);X(num_X,6)], NAV(i).rxTime);
    else
        relp = [NAV(1).posX;NAV(1).posY;NAV(1).posZ];
        velp = [NAV(1).velX;NAV(1).velY;NAV(1).velZ];
    end

    % [pos(i,:), dt(i), idx, dpos ] = DPE_rinex_Moon(satPos, Grid_ECEF, Origin_ECEF);
    [pos(i,:), dt(i), idx, dpos ] = DPE_rinex_Moon_orbit3(satPos, Grid_ECEF, Origin_ECEF, relp, i);
    % [pos(i,:), dt(i), idx, dpos ] = DPE_rinex_Moon_orbit_vel(satPos, Grid_ECEF, Origin_ECEF, velp, dtime, i);

    % [pos2(i,:), dt2(i), dop, others] = leastSquarePos(satPos);

    % [pos3(i,:), dt3(i), idx2, dpos2] = DPE_rinex_Moon_dop(satVel, NAV(i).velX, NAV(i).velY, NAV(i).velZ, Grid_ECEF, Origin_ECEF);

    dx(i) = abs(pos(i,1) - NAV(i).posX);
    dy(i) = abs(pos(i,2) - NAV(i).posY);
    dz(i) = abs(pos(i,3) - NAV(i).posZ);
    dd(i) = sqrt(dx(i)^2+dy(i)^2+dz(i)^2);

    % dx(i) = abs(pos(i,1) - pos2(i,1));
    % dy(i) = abs(pos(i,2) - pos2(i,2));
    % dz(i) = abs(pos(i,3) - pos2(i,3));
    % dd(i) = sqrt(dx(i)^2+dy(i)^2+dz(i)^2);

    dx2(i) = abs(relp(1) - NAV(i).posX);
    dy2(i) = abs(relp(2) - NAV(i).posY);
    dz2(i) = abs(relp(3) - NAV(i).posZ);
    dd2(i) = sqrt(dx2(i)^2+dy2(i)^2+dz2(i)^2);

    dx3(i) = abs(pos2(i,1) - NAV(i).posX);
    dy3(i) = abs(pos2(i,2) - NAV(i).posY);
    dz3(i) = abs(pos2(i,3) - NAV(i).posZ);
    dd3(i) = sqrt(dx3(i)^2+dy3(i)^2+dz3(i)^2);

    % dx3(i) = abs(pos3(i,1) - NAV(i).posX);
    % dy3(i) = abs(pos3(i,2) - NAV(i).posY);
    % dz3(i) = abs(pos3(i,3) - NAV(i).posZ);
    % dd3(i) = sqrt(dx3(i)^2+dy3(i)^2+dz3(i)^2);

    Origin_ECEF = pos(i,:)';
    % if i < n
    %     Origin_ECEF = [NAV(i+1).posX; NAV(i+1).posY; NAV(i+1).posZ];
    % end
    dt0 = dt(i);
    [Xr, Vr] = ecef_to_j2000(pos(count_num2,:)', [NAV(count_num2).velX;NAV(count_num2).velY;NAV(count_num2).velZ], NAV(count_num2).rxTime);
    X0 = [Xr;Vr;moonEph(floor((i-1)/60)+1,2);moonEph(floor((i-1)/60)+1,3);moonEph(floor((i-1)/60)+1,4);moonEph(floor((i-1)/60)+1,5);moonEph(floor((i-1)/60)+1,6);moonEph(floor((i-1)/60)+1,7)];
    % [Xr, Vr] = ecef_to_j2000([NAV(count_num2).posX;NAV(count_num2).posY;NAV(count_num2).posZ], [NAV(count_num2).velX;NAV(count_num2).velY;NAV(count_num2).velZ], NAV(count_num2).rxTime);
    % X0 = [Xr;Vr;moonEph(floor((i-1)/60)+1,2);moonEph(floor((i-1)/60)+1,3);moonEph(floor((i-1)/60)+1,4);moonEph(floor((i-1)/60)+1,5);moonEph(floor((i-1)/60)+1,6);moonEph(floor((i-1)/60)+1,7)];
    str=['Processing ... ',num2str(i), ' / ', num2str(n)];
    waitbar(i/n, hwait, str);
end
close(hwait);
fprintf('\n');


%% ===== Percentile statistics (Position error) =====

% pList = [25 50 75 95];

% 去掉 NaN（非常重要）
dd_valid  = dd(~isnan(dd));
% dd2_valid = dd2(~isnan(dd2));
dd3_valid = dd3(~isnan(dd3));

% % DPE percentiles
% DPE_percentiles = prctile(dd_valid, pList);
% 
% % Least Squares percentiles
% LS_percentiles  = prctile(dd2_valid, pList);
% 
% fprintf('\n===== Position Error Percentiles (meters) =====\n');
% 
% fprintf('DPE Method:\n');
% fprintf('  25%% : %.3f m\n', DPE_percentiles(1));
% fprintf('  50%% : %.3f m (median)\n', DPE_percentiles(2));
% fprintf('  75%% : %.3f m\n', DPE_percentiles(3));
% fprintf('  95%% : %.3f m\n', DPE_percentiles(4));
% 
% fprintf('\nLeast Squares Method:\n');
% fprintf('  25%% : %.3f m\n', LS_percentiles(1));
% fprintf('  50%% : %.3f m (median)\n', LS_percentiles(2));
% fprintf('  75%% : %.3f m\n', LS_percentiles(3));
% fprintf('  95%% : %.3f m\n', LS_percentiles(4));



diffpos = diff(pos,1,1);
diffdt = diff(dt,1);

mindt = min(diffdt,[],'omitnan');
maxdt = max(diffdt,[],'omitnan');

DPExleft  = min(diffpos(:,1),[],'omitnan');
DPExright = max(diffpos(:,1),[],'omitnan');
DPEyleft  = min(diffpos(:,2),[],'omitnan');
DPEyright = max(diffpos(:,2),[],'omitnan');
DPEzleft  = min(diffpos(:,3),[],'omitnan');
DPEzright = max(diffpos(:,3),[],'omitnan');

disp(['DPEx: ',  num2str(DPExleft), ' ', num2str(DPExright)]);
disp(['DPEy: ', num2str(DPEyleft), ' ', num2str(DPEyright)]);
disp(['DPEz: ', num2str(DPEzleft), ' ', num2str(DPEzright)]);
disp(['DPEdt: ', num2str(mindt), ' ', num2str(maxdt)]);

figure(1)
plot(1:n,dx)
xlabel('dx (m)');
ylabel('error (m)');

figure(2)
plot(1:n,dy)
xlabel('dy (m)');
ylabel('error (m)');

figure(3)
plot(1:n,dz)
xlabel('dz (m)');
ylabel('error (m)');

figure(4)
plot(1:n,dt)
ylabel('钟差/m')

figure(5)
plot(1:n,dd)
ylabel('定位误差/m')

figure(6)
plot(1:n,dd2)
ylabel('定位误差/m')

figure(7)
plot(1:n,dd3)
ylabel('定位误差/m')

figure(8)  % 使用同一个图窗编号
plot(1:n, dd, 'b-', 'LineWidth', 1.5);  % dd 用蓝色实线
hold on;
plot(1:n, dd3, 'r-', 'LineWidth', 1.5);  % dd3 用红色实线
hold off;

xlabel('历元');
ylabel('3D定位误差 / m');
legend('dd', 'dd3');  % 添加图例
% grid on;  % 可选：添加网格，便于观察
% disp(RAW.signalId(1))

% ===== 画 3D 散点图 =====
X2 = zeros(n,1);
Y2 = zeros(n,1);
Z2 = zeros(n,1);

for i = 1:n
    X2(i) = NAV(i).posX;
    Y2(i) = NAV(i).posY;
    Z2(i) = NAV(i).posZ;
end

figure;
hold on;
grid on;
axis equal;

% 第一组：n×3 数组（红色）
scatter3(pos(:,1), pos(:,2), pos(:,3), ...
         40, 'r', 'filled');

scatter3(pos2(:,1), pos2(:,2), pos2(:,3), ...
         40, 'g', 'filled');

% 第二组：struct 中的坐标（蓝色）
scatter3(X2, Y2, Z2, ...
         40, 'b', 'filled');

xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
title('ECEF Coordinates Comparison');

legend('DPE Results', 'LS Results' ,'References');

view(3);