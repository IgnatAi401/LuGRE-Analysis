function [gpsData, galileoData] = read_gnss_rinex(filename)
    % 1. 打开文件
    fid = fopen(filename, 'r');
    if fid == -1
        error('无法打开文件: %s', filename);
    end

    % ---------------------------------------------------------
    % 预估最大历元数量并预分配内存
    % ---------------------------------------------------------
    MAX_GUESS = 5000; % 初始预估大小
    
    % 定义空的结构体模板
    % PRN 初始化为 0 (Double类型)
    dummyStruct = struct(...
        'PRN', 0, 'Time', datetime, 'GPSTotalSeconds', 0, ...
        'ClockBias', 0, 'ClockDrift', 0, 'ClockDriftRate', 0, ...
        'IODE', 0, 'Crs', 0, 'DeltaN', 0, 'M0', 0, ...
        'Cuc', 0, 'e', 0, 'Cus', 0, 'SqrtA', 0, ...
        'Toe', 0, 'Cic', 0, 'Omega0', 0, 'Cis', 0, ...
        'i0', 0, 'Crc', 0, 'omega', 0, 'OmegaDot', 0, ...
        'IDOT', 0, 'CodesL2', 0, 'Week', 0, 'L2PFlag', 0, ...
        'SVAccuracy', 0, 'SVHealth', 0, 'TGD', 0, 'IODC', 0, ...
        'TxTime', 0);

    % 预分配结构体数组
    gpsTemp = repmat(dummyStruct, MAX_GUESS, 1);
    galTemp = repmat(dummyStruct, MAX_GUESS, 1);
    
    g_count = 0; 
    e_count = 0; 

    % 计算 GPS 起始时间 (1980-01-06) 用于计算总秒数
    GPS_EPOCH_NUM = datenum(1980, 1, 6, 0, 0, 0);

    % 4. 跳过 Header
    while ~feof(fid)
        line = fgetl(fid);
        if contains(line, 'END OF HEADER')
            break;
        end
    end

    % 5. 循环读取数据体
    while ~feof(fid)
        line = fgetl(fid);
        if ~ischar(line) || length(line) < 2
            continue;
        end
        
        sysID = line(1);
        
        switch sysID
            case 'G' % GPS
                g_count = g_count + 1;
                % 如果超出预分配大小，自动倍增扩展
                if g_count > length(gpsTemp)
                    gpsTemp = [gpsTemp; repmat(dummyStruct, MAX_GUESS, 1)]; 
                end
                block = read_block_fast(fid, line, 7);
                gpsTemp(g_count) = parse_nav_record_fast(block, GPS_EPOCH_NUM);
                
            case 'E' % Galileo
                e_count = e_count + 1;
                if e_count > length(galTemp)
                    galTemp = [galTemp; repmat(dummyStruct, MAX_GUESS, 1)];
                end
                block = read_block_fast(fid, line, 7);
                galTemp(e_count) = parse_nav_record_fast(block, GPS_EPOCH_NUM);
                
            case {'C', 'J', 'I'} 
                read_block_fast(fid, line, 7); 
            case {'R', 'S'}
                read_block_fast(fid, line, 3); 
        end
    end
    fclose(fid);

    % 6. 裁剪结果 (返回给函数输出变量)
    gpsData = gpsTemp(1:g_count);
    galileoData = galTemp(1:e_count);

    % 7. 简单的结果展示 (可选)
    if ~isempty(gpsData)
        disp('读取成功。');
        disp(['GPS 记录数: ', num2str(length(gpsData))]);
        disp(['Galileo 记录数: ', num2str(length(galileoData))]);
        % 打印第一条示例验证 PRN 和 时间
        disp(['示例 GPS PRN: ', num2str(gpsData(1).PRN), ...
              ' (TotalSec: ', num2str(gpsData(1).GPSTotalSeconds, '%.1f'), ')']);
    end
end

%% 辅助函数：快速读取块
function block = read_block_fast(fid, firstLine, numDataLines)
    block = cell(1 + numDataLines, 1);
    block{1} = firstLine;
    for i = 1:numDataLines
        block{1+i} = fgetl(fid);
    end
end

%% 辅助函数：快速解析
function rec = parse_nav_record_fast(lines, gps_epoch_num)
    line1 = lines{1};
    
    % 解析 PRN 为 double
    % RINEX 格式取后两位转数字
    rec.PRN = str2double(line1(2:3));
    
    % 时间解析
    year = str2double(line1(5:8));
    month = str2double(line1(10:11));
    day = str2double(line1(13:14));
    hour = str2double(line1(16:17));
    minute = str2double(line1(19:20));
    second = str2double(line1(22:23));
    
    % 保存 datetime
    rec.Time = datetime(year, month, day, hour, minute, second);
    
    % [核心需求] 计算 GPS 总秒数
    current_datenum = datenum(year, month, day, hour, minute, second);
    rec.GPSTotalSeconds = (current_datenum - gps_epoch_num) * 86400;
    
    % 读取钟差
    rec.ClockBias = str2num_fast(line1(24:42));
    rec.ClockDrift = str2num_fast(line1(43:61));
    rec.ClockDriftRate = str2num_fast(line1(62:80));
    
    % 展开读取剩余7行参数
    l2 = pad(lines{2}, 80);
    rec.IODE = str2num_fast(l2(5:23));
    rec.Crs = str2num_fast(l2(24:42));
    rec.DeltaN = str2num_fast(l2(43:61));
    rec.M0 = str2num_fast(l2(62:80));
    
    l3 = pad(lines{3}, 80);
    rec.Cuc = str2num_fast(l3(5:23));
    rec.e = str2num_fast(l3(24:42));
    rec.Cus = str2num_fast(l3(43:61));
    rec.SqrtA = str2num_fast(l3(62:80));
    
    l4 = pad(lines{4}, 80);
    rec.Toe = str2num_fast(l4(5:23));
    rec.Cic = str2num_fast(l4(24:42));
    rec.Omega0 = str2num_fast(l4(43:61));
    rec.Cis = str2num_fast(l4(62:80));
    
    l5 = pad(lines{5}, 80);
    rec.i0 = str2num_fast(l5(5:23));
    rec.Crc = str2num_fast(l5(24:42));
    rec.omega = str2num_fast(l5(43:61));
    rec.OmegaDot = str2num_fast(l5(62:80));
    
    l6 = pad(lines{6}, 80);
    rec.IDOT = str2num_fast(l6(5:23));
    rec.CodesL2 = str2num_fast(l6(24:42));
    rec.Week = str2num_fast(l6(43:61));
    rec.L2PFlag = str2num_fast(l6(62:80));
    
    l7 = pad(lines{7}, 80);
    rec.SVAccuracy = str2num_fast(l7(5:23));
    rec.SVHealth = str2num_fast(l7(24:42));
    rec.TGD = str2num_fast(l7(43:61));
    rec.IODC = str2num_fast(l7(62:80));
    
    l8 = pad(lines{8}, 80);
    rec.TxTime = str2num_fast(l8(5:23));
end

%% 辅助函数：安全的字符串转数字
function val = str2num_fast(str)
    if isempty(str)
        val = 0; return;
    end
    if contains(str, 'D')
        str = strrep(str, 'D', 'E');
    end
    val = str2double(str);
    if isnan(val), val = 0; end
end