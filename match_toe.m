function RAW_sorted = match_toe(RAW_sorted, csv_filepath)
% MATCH_TOE_DIRECT 匹配CSV中的星历数据并将最近的Toe直接赋值给结构体
%
% 输入:
%   RAW_sorted   - 1xn 的结构体，包含 txTime, signalId, svId
%   csv_filepath - .csv 文件的路径字符串
%
% 输出:
%   RAW_sorted   - 更新后的结构体，新增字段 Toe (直接取自CSV，未做修正)
%

    %% 1. 读取 CSV 数据
    fprintf('正在读取 CSV 文件...\n');
    try
        % 假设数据没有标题行，或者readmatrix能自动识别数值区域
        data = readmatrix(csv_filepath);
    catch ME
        error('读取CSV文件失败: %s', ME.message);
    end
    
    % 提取 CSV 各列数据
    csv_rxTime = data(:, 1);
    csv_signal = data(:, 2);
    csv_svId   = data(:, 3);
    csv_Toe    = data(:, 4); % 我们只需要这一列的值
    
    %% 2. 准备结构体数据 (向量化预处理)
    fprintf('正在匹配数据...\n');
    
    % 将结构体字段转换为向量，提高索引速度
    s_txTime = [RAW_sorted.txTime];
    s_signalId = [RAW_sorted.signalId];
    s_svId = [RAW_sorted.svId];
    
    % 初始化结果向量（默认填充 NaN，防止找不到匹配时报错）
    n_samples = length(s_txTime);
    s_Toe_new = nan(1, n_samples);
    
    %% 3. 按卫星 (Signal + SvId) 分组处理
    % 找出 RAW_sorted 中出现过的所有唯一的 [signal, svId] 组合
    unique_pairs = unique([s_signalId', s_svId'], 'rows');
    
    for i = 1:size(unique_pairs, 1)
        % 当前处理的卫星标识
        u_sig = unique_pairs(i, 1);
        u_sv  = unique_pairs(i, 2);
        
        % A. 筛选 CSV 数据：找到该卫星的所有行
        idx_csv = (csv_signal == u_sig) & (csv_svId == u_sv);
        if ~any(idx_csv)
            continue; % 如果CSV里没有这颗卫星的数据，跳过
        end
        
        % 提取该卫星在 CSV 中的时间 和 Toe
        curr_csv_rxTime = csv_rxTime(idx_csv);
        curr_csv_Toe    = csv_Toe(idx_csv);
        
        % B. 筛选结构体数据：找到该卫星的所有待处理数据
        idx_struct = (s_signalId == u_sig) & (s_svId == u_sv);
        
        % 提取结构体中待处理的 txTime
        curr_struct_txTime = s_txTime(idx_struct);
        
        % C. 寻找最近的时间点
        % dsearchn 用于在 curr_csv_rxTime 中找到离 curr_struct_txTime 最近的点的索引
        % nearest_indices 对应的是 curr_csv_Toe 中的下标
        nearest_indices = dsearchn(curr_csv_rxTime, curr_struct_txTime(:));
        
        % D. 直接赋值
        % 找到索引对应的 Toe 值，直接填入结果向量
        s_Toe_new(idx_struct) = curr_csv_Toe(nearest_indices);
    end
    
    %% 4. 将结果写回结构体
    % 使用 num2cell 和 deal 快速将向量分配给结构体字段
    toe_cell = num2cell(s_Toe_new);
    [RAW_sorted.Toe] = toe_cell{:};
    
    fprintf('处理完成。\n');
end