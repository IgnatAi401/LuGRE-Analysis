clc; clear;
basePath = '/Users/1gnat4y/Downloads/LuGRE';

OP_tags = {
    'OP1_0'
    'OP2_0'
    'OP5_0'
    'OP9_0'
    'OP14_0'
    'OP17_0'
    'OP21_0'
    'OP38_0'
    'OP40_0'
    'OP74_0'
    'OP76_0'
    'OP77_0'
    'OP77_1'
    'OP78_1'
};

ACQ_all = struct([]);
RAW_all = struct([]);
NAV_all = struct([]);

for k = 1:numel(OP_tags)
    OP_str = OP_tags{k};

    % ========= 新增：根据 OP_tag 判定 type =========
    if ismember(OP_str, {'OP1_0', 'OP2_0'})
        type_val = 0;
    elseif ismember(OP_str, {'OP5_0','OP9_0','OP14_0','OP17_0','OP21_0'})
        type_val = 1;
    else
        % OP38 ~ OP78 统一设为 3
        type_val = 2;
    end
    % ===============================================

    [ACQfilePath, RAWfilePath, NAVfilePath] = ...
        FindTelemetryTxtFiles(basePath, OP_str);

    ACQ = TxtParserAcq(ACQfilePath);
    RAW = TxtParserRaw(RAWfilePath);
    NAV = TxtParserNav(NAVfilePath);

    % 标注来源 + type
    for i = 1:numel(ACQ)
        ACQ(i).OP_tag = OP_str;
        ACQ(i).type   = type_val;
    end
    for i = 1:numel(RAW)
        RAW(i).OP_tag = OP_str;
        RAW(i).type   = type_val;
    end
    for i = 1:numel(NAV)
        NAV(i).OP_tag = OP_str;
        NAV(i).type   = type_val;
    end

    ACQ_all = [ACQ_all, ACQ];
    RAW_all = [RAW_all, RAW];
    NAV_all = [NAV_all, NAV];
end
