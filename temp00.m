% 假设 RAW 是 1*n 的结构体，包含字段 rxTime, signalId, svId, cn0

%% 1. 输出所有不重复的 (signalId, svId) 组合
comb = unique([ [RAW.signalId]' , [RAW.svId]' ], 'rows');
disp('所有不重复的 signalId + svId 组合：');
disp(array2table(comb, 'VariableNames', {'signalId', 'svId'}));

% （可选）如果你想看到信号名称，可以添加映射：
signalNames = {'GPS_L1CA', 'GPS_L5', 'GAL_E1BC', 'GAL_E5A', 'GAL_E5B'};
fprintf('\n组合详情（含信号名）：\n');
for i = 1:size(comb, 1)
    fprintf('  signalId=%d (%s), svId=%d\n', ...
        comb(i,1), signalNames{comb(i,1)+1}, comb(i,2));
end

%% 2. 统计 cn0 的分布
cn0 = [RAW.cn0];
fprintf('\n===== CN0 统计分布 =====\n');
fprintf('最小值     : %.2f dB-Hz\n', min(cn0));
fprintf('最大值     : %.2f dB-Hz\n', max(cn0));
fprintf('平均值     : %.2f dB-Hz\n', mean(cn0));
fprintf('中位数     : %.2f dB-Hz\n', median(cn0));
fprintf('标准差     : %.2f dB-Hz\n', std(cn0));
p = prctile(cn0, [25, 75]);
fprintf('25分位数   : %.2f dB-Hz\n', p(1));
fprintf('75分位数   : %.2f dB-Hz\n', p(2));
fprintf('总观测量数 : %d\n', length(cn0));

% 画直方图
figure;
histogram(cn0, 'BinWidth', 1);  % 可按需要调整 BinWidth
xlabel('C/N_0 (dB-Hz)');
ylabel('频数');
title('cn0 分布直方图');
grid on;

%% 3. 将 rxTime 的开始和结束时间转换为人类可读的格式
rxAll = [RAW.rxTime];
gpsEpoch = datetime(1980, 1, 6, 0, 0, 0); % GPS 起始历元 (UTC)

% 注意：GPS 时间比 UTC 快若干闰秒（2026年时快 18 秒）。
% 这里简单以 GPS 时间直接转为日历时间，若需要精确 UTC 可减去 18 秒：
% startDateTime = gpsEpoch + seconds(min(rxAll)) - seconds(18);
startDateTime = gpsEpoch + seconds(min(rxAll));
endDateTime   = gpsEpoch + seconds(max(rxAll));

fprintf('\n===== 观测时间范围 =====\n');
fprintf('开始时间: %s\n', datestr(startDateTime, 'yyyy-mm-dd HH:MM:SS'));
fprintf('结束时间: %s\n', datestr(endDateTime,   'yyyy-mm-dd HH:MM:SS'));