% clear; clc;

% t_gps_array = zeros(1,1387);
% for i = 1:1387
%     t_gps_array(i) = NAV(i).rxTime;
% end



function results = getMoonEph(t_gps_array)
    % ===== 1. 时间准备 =====
    gps_epoch = datetime(1980,1,6,0,0,0);
    t_start = gps_epoch + seconds(min(t_gps_array));
    t_end   = gps_epoch + seconds(max(t_gps_array));
    
    
    % 关键：NASA 接受 'YYYY-MM-DD HH:MM:SS' 格式
    t_start_str = datestr(t_start, 'yyyy-mm-dd HH:MM:SS');
    t_end_str   = datestr(t_end,   'yyyy-mm-dd HH:MM:SS');
    
    % 
    %     % ===== 2. API请求 =====
    %     baseUrl = 'https://ssd.jpl.nasa.gov/api/horizons.api';
    % 
    %     % 注意：这里去掉了所有多余的 ''' 嵌套引号
    %     % webread 会自动处理空格和特殊字符的转义 (URL Encoding)
    %     params = struct(...
    %         'format', 'text', ...
    %         'COMMAND', '301', ...    
    %         'OBJ_DATA', 'NO', ...    
    %         'MAKE_EPHEM', 'YES', ... 
    %         'EPHEM_TYPE', 'VECTORS', ...
    %         'CENTER', '500@0', ...   
    %         'REF_SYSTEM', 'J2000', ...
    %         'OUT_UNITS', 'KM-S', ...
    %         'TIME_TYPE', 'UTC', ...
    %         'VEC_TABLE', '2', ...
    %         'START_TIME', t_start_str, ... 
    %         'STOP_TIME', t_end_str, ...
    %         'STEP_SIZE', step_str);
    % 
    %     options = weboptions('Timeout', 60); % 增加超时时间
    % 
    %     try
    %         data = webread(baseUrl, params, options);
    %         disp('数据请求成功！');
    %     catch ME
    %         warning('API请求失败！错误信息：%s', ME.message);
    %         % 如果失败了，在这里提前退出，返回空值
    %         return; 
    %     end
    % 
    %     % ===== 3. 解析数据 =====
    %     lines = splitlines(string(data));
    %     idx_start = find(contains(lines, '$$SOE')) + 1;
    %     idx_end   = find(contains(lines, '$$EOE')) - 1;
    % 
    %     if isempty(idx_start) || isempty(idx_end)
    %         warning('未能找到历表数据标识符 $$SOE/$$EOE。请检查 API 返回内容。');
    %         disp(data); % 打印出来看 API 返回了什么错误
    %         return;
    %     end
    % 
    %     % 提取坐标和速度
    %     % Horizons VEC_TABLE=2 的格式通常是每组 3 行
    %     n_data = floor((idx_end - idx_start + 1) / 3);
    %     r_all = zeros(n_data, 3);
    %     v_all = zeros(n_data, 3);
    % 
    %     for k = 1:n_data
    %         base_line = idx_start + (k-1)*3;
    %         % 解析位置 X, Y, Z
    %         tokens_r = regexp(lines(base_line+1), 'X =\s*([-\d.E+]+)\s*Y =\s*([-\d.E+]+)\s*Z =\s*([-\d.E+]+)', 'tokens');
    %         if ~isempty(tokens_r)
    %             r_all(k,:) = str2double(tokens_r{1});
    %         end
    %         % 解析速度 VX, VY, VZ
    %         tokens_v = regexp(lines(base_line+2), 'VX=\s*([-\d.E+]+)\s*VY=\s*([-\d.E+]+)\s*VZ=\s*([-\d.E+]+)', 'tokens');
    %         if ~isempty(tokens_v)
    %             v_all(k,:) = str2double(tokens_v{1});
    %         end
    %     end
    % 
    %     % ===== 4. 插值到原时间 =====
    %     % 生成 API 返回的对应时间序列
    %     t_api_points = linspace(min(t_gps_array), max(t_gps_array), n_data);
    %     r_moon = interp1(t_api_points, r_all, t_gps_array, 'spline');
    %     v_moon = interp1(t_api_points, v_all, t_gps_array, 'spline');
    % end
    
    
    
    
    %% 参数设置
    target = '301';           
    center = '500@399';       
    % 注意：NASA 建议时间格式为 'YYYY-MMM-DD HH:MM:SS'，例如 '2026-Mar-24 12:00:00'
    % time_point = '2026-03-21 12:00:00'; 
    % 
    % time_point2 = '2026-03-21 15:00:00'; 
    
    base_url = 'https://ssd.jpl.nasa.gov/api/horizons.api';
    query_url = [base_url, ...
        '?format=json', ...
        '&COMMAND=''', target, '''', ...
        '&OBJ_DATA=''NO''', ...
        '&MAKE_EPHEM=''YES''', ...
        '&EPHEM_TYPE=''VECTORS''', ...
        '&CENTER=''', center, '''', ...
        '&START_TIME=''', t_start_str, '''', ...
        '&STOP_TIME=''', t_end_str, '''', ... % 多给1分钟确保覆盖
        '&STEP_SIZE=''1m''', ... 
        '&VEC_TABLE=''2''', ...
        '&REF_PLANE=''ECLIPTIC''', ... 
        '&OUT_UNITS=''KM-S'''];
    
    %% 发送请求
    
    options = weboptions('Timeout', 60);
    response = webread(query_url, options);
    raw_data = response.result;
    
    % 检查是否存在数据标记
    soe_idx = strfind(raw_data, '$$SOE');
    eoe_idx = strfind(raw_data, '$$EOE');
    
    if isempty(soe_idx) || isempty(eoe_idx)
        fprintf('--- API 返回错误信息 ---\n');
        disp(raw_data); % 打印 NASA 返回的具体报错原因
        error('未能从返回结果中定位到数据段 ($$SOE)。');
    end
    
    % 提取第一行数据
    content = raw_data(soe_idx+5 : eoe_idx-1);
    % 假设 a 是原始长字符串
    
    % 1. 预分配空间（根据你的数据量估算，或者动态增长）
    % 我们利用正则表达式把数据切分成每一个时间点的"小块"
    % 每一个块都以 Julian Date 开头
    blocks = regexp(content, '\d{7}\.\d{9} = A\.D\..*?(?=\d{7}\.\d{9} = A\.D\.|$)', 'match');
    
    n = length(blocks);
    results = zeros(n, 7); % 存储 [JD, X, Y, Z, VX, VY, VZ]
    
    for i = 1:n
        current_block = blocks{i};
        
        % 提取 JD
        jd_val = regexp(current_block, '^(\d+\.\d+)', 'tokens', 'once');
        if ~isempty(jd_val), results(i, 1) = str2double(jd_val{1}); end
        
        % 提取 X, Y, Z (精确匹配，确保后面紧跟空格或换行)
        x_val = regexp(current_block, 'X\s*=\s*([\d\.E\+\-]+)', 'tokens', 'once');
        y_val = regexp(current_block, 'Y\s*=\s*([\d\.E\+\-]+)', 'tokens', 'once');
        z_val = regexp(current_block, 'Z\s*=\s*([\d\.E\+\-]+)', 'tokens', 'once');
        
        % 提取 VX, VY, VZ (明确指定匹配 VX)
        vx_val = regexp(current_block, 'VX\s*=\s*([\d\.E\+\-]+)', 'tokens', 'once');
        vy_val = regexp(current_block, 'VY\s*=\s*([\d\.E\+\-]+)', 'tokens', 'once');
        vz_val = regexp(current_block, 'VZ\s*=\s*([\d\.E\+\-]+)', 'tokens', 'once');
        
        % 填充数据
        if ~isempty(x_val),  results(i, 2) = str2double(x_val{1});  end
        if ~isempty(y_val),  results(i, 3) = str2double(y_val{1});  end
        if ~isempty(z_val),  results(i, 4) = str2double(z_val{1});  end
        if ~isempty(vx_val), results(i, 5) = str2double(vx_val{1}); end
        if ~isempty(vy_val), results(i, 6) = str2double(vy_val{1}); end
        if ~isempty(vz_val), results(i, 7) = str2double(vz_val{1}); end
    end
    
    % % 转换为 Table 方便观察
    % final_table = array2table(results, 'VariableNames', ...
    %     {'JD', 'X_km', 'Y_km', 'Z_km', 'VX_kms', 'VY_kms', 'VZ_kms'});
    
    % % 导出 Excel
    % writetable(final_table, 'Processed_Ephemeris.xlsx');
    % disp('数据已清洗完成，无重复穿插。');
end