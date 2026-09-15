%% 1. 数据检查与提取
clc; clear;
load total.mat

fprintf('正在处理 %d 条观测数据...\n', length(RAW_all));

% 向量化提取字段
all_rxTime   = [RAW_all.rxTime];
all_signalId = [RAW_all.signalId];
all_svId     = [RAW_all.svId];
all_cn0      = [RAW_all.cn0];

% 相对时间（秒）
time_plot = all_rxTime - all_rxTime(1); 

%% 2. 定义信号类型映射
% 0:GPS_L1CA, 1:GPS_L5, 2:GAL_E1BC, 3:GAL_E5A, 4:GAL_E5B
sigMap = containers.Map({0, 1, 2, 3, 4}, ...
    {'GPS L1CA', 'GPS L5', 'GAL E1BC', 'GAL E5A', 'GAL E5B'});

unique_signals = unique(all_signalId);

%% 3. 开始绘图（散点图）
figure('Name', 'GNSS C/N0 Scatter Plot', ...
       'Color', 'w', 'NumberTitle', 'off');
tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:length(unique_signals)
    curr_sig = unique_signals(i);
    
    % --- A: 当前信号类型 ---
    idx_sig = (all_signalId == curr_sig);
    
    t_sub   = time_plot(idx_sig);
    sv_sub  = all_svId(idx_sig);
    cn0_sub = all_cn0(idx_sig);
    
    unique_svs = unique(sv_sub);
    
    % 子图
    nexttile;
    hold on; grid on; box on;
    
    colors = lines(length(unique_svs));
    
    for k = 1:length(unique_svs)
        prn = unique_svs(k);
        idx_prn = (sv_sub == prn);
        
        t_prn = t_sub(idx_prn);
        c_prn = cn0_sub(idx_prn);
        
        % ===== 核心修改：散点图 =====
        scatter(t_prn, c_prn, ...
            18, ...                         % 点大小
            colors(k,:), ...               % 颜色
            'filled', ...
            'DisplayName', sprintf('PRN %d', prn));
    end
    
    % --- B: 图表美化 ---
    if isKey(sigMap, curr_sig)
        titleStr = sigMap(curr_sig);
    else
        titleStr = sprintf('Unknown Signal ID %d', curr_sig);
    end
    
    title(titleStr, 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('Relative Time (s)');
    ylabel('C/N0 (dB-Hz)');
    
    ylim([20 60]);
    
    legend('show', ...
           'Location', 'bestoutside', ...
           'Interpreter', 'none');
end

sgtitle('C/N0 Scatter Distribution by Signal Type');
