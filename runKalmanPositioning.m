% 运行前要求工作区中已经存在 NAV 和 RAW_Pos。
% getPos.m 也需要位于 MATLAB 当前路径中。

if ~exist('NAV', 'var') || ~exist('RAW_Pos', 'var')
    error('工作区中必须先存在 NAV 和 RAW_Pos 变量。');
end

cfg = struct();
cfg.dt = 1.0;

% 默认采用 GPS L1。如果第 9 列不是 GPS L1 Doppler，请改为实际载波频率。
cfg.carrierFrequencyHz = 1575.42e6;

% 常见定义：伪距率 = -lambda * Doppler。
% 如果滤出的速度明显反向或异常，可尝试改为 +1。
cfg.dopplerSign = -1;

% 静态接收机可尝试 0.01~0.1；车辆一般可用 0.5~5。
cfg.accelPSD = 0.5;

kfResult = kalmanGnssPositioning(NAV, RAW_Pos, cfg);

valid = kfResult.initialized;
fprintf('完成：共 %d 个历元，成功输出 %d 个历元。\n', ...
    numel(NAV), sum(valid));

if any(valid)
    figure('Name', 'GNSS EKF Position');
    subplot(3, 1, 1);
    plot(kfResult.posECEF(:, 1)); grid on;
    ylabel('X / m');
    subplot(3, 1, 2);
    plot(kfResult.posECEF(:, 2)); grid on;
    ylabel('Y / m');
    subplot(3, 1, 3);
    plot(kfResult.posECEF(:, 3)); grid on;
    ylabel('Z / m'); xlabel('Epoch');

    figure('Name', 'GNSS EKF Velocity');
    plot(kfResult.velECEF, 'LineWidth', 1.0); grid on;
    xlabel('Epoch'); ylabel('Velocity / (m/s)');
    legend('V_x', 'V_y', 'V_z', 'Location', 'best');

    finalEpoch = find(valid, 1, 'last');
    fprintf('末历元 ECEF 位置 (m): %.3f, %.3f, %.3f\n', ...
        kfResult.posECEF(finalEpoch, :));
    fprintf('末历元 LLH: lat %.9f deg, lon %.9f deg, h %.3f m\n', ...
        kfResult.posLLH(finalEpoch, :));
end
