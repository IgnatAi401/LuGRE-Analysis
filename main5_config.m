function cfg = main5_config()
% MAIN5_CONFIG main5 的统一配置。距离 m，速度 m/s，钟差 c*dt (m)。
root = fileparts(mfilename('fullpath'));

%% 数据与运行范围
cfg.data.root = '/Users/ignat/Documents/LuGRE';
cfg.data.window = 'OP1_0';
cfg.data.rinex = fullfile(root, 'BRDC00IGS_R_20250150000_01D_MN.rnx');
cfg.run.maxEpochs = Inf;             % 调试时可设 5、60；Inf 表示全部
cfg.run.plot = true;
cfg.run.progressEvery = 60;         % 每多少历元打印进度，0 表示关闭

%% 方法开关（true 开启，false 关闭）
cfg.enable.dpe = true;
cfg.enable.ls = true;
cfg.enable.cn0Dpe = false;           % DPE 残差按 C/N0 加权
cfg.enable.cn0Ls = false;            % LS 残差平方按 C/N0 加权
cfg.enable.orbit = true;            % 轨道传播 + DPE 位置先验；LS 仍为纯观测对照
cfg.enable.j2 = true;               % 仅 orbit=true 时生效
cfg.enable.sagnac = true;           % 地球自转传播修正，两种方法共用
cfg.enable.legacyIonosphere = false;% 原经验公式不适合地月几何，仅供复现旧结果

%% 初始状态（继续使用 NAV 初始位置、分段 NAV 速度）
cfg.initial.positionOffsetM = [0; 0; 0];
cfg.initial.clockBiasM = -140000;

%% DPE 搜索网格：相对上一有效解的偏移
cfg.dpe.switchEpoch = 201;
cfg.dpe.coarseOffsetsM = -20000:1000:20000;
cfg.dpe.fineOffsetsM = -10000:500:10000;
cfg.dpe.clockOffsetsM = -25000:1000:10000;
cfg.dpe.orbitWeight = 1e-5;          % 保留原值；距离惩罚系数，尚未标定
cfg.dpe.chunkSize = 200000;         % 分块计算候选点，降低内存峰值

%% LS 迭代
cfg.ls.maxIterations = 30;
cfg.ls.positionToleranceM = 0.01;
cfg.ls.clockToleranceM = 0.01;
cfg.ls.rankTolerance = 1e-10;

%% 观测筛选与信噪比权重
cfg.observation.minCount = 4;
cfg.observation.cn0ReferenceDbHz = 30;
cfg.observation.weightLimits = [0.01, 100]; % w=10^((C/N0-reference)/10)，再限幅
cfg.observation.legacyZenithDelayS = 5e-9;

%% 轨道传播（地球 + 月球 + 可选地球 J2）
cfg.orbit.resetIntervalS = 60;      % 按实际时间重置，不按样本编号
cfg.orbit.relativeTolerance = 1e-9;
cfg.orbit.positionAbsToleranceM = 0.01;
cfg.orbit.velocityAbsToleranceMps = 1e-5;

%% 时间、坐标系与月球星历
cfg.frame.gpsMinusUtcS = 18;        % 本项目 2025 年数据；换年代需核对闰秒
cfg.frame.eopSource = 'iers';       % 'iers': MATLAB IERS 数据；'manual': 使用下两项
cfg.frame.deltaUt1S = 0;            % manual 模式：UT1-UTC (s)
cfg.frame.polarMotionRad = [0, 0];  % manual 模式：[xp yp] (rad)，0 是近似
cfg.frame.derivativeStepS = 0.5;    % 旋转矩阵时间导数，供速度转换
cfg.moon.stepMinutes = 1;
cfg.moon.timeoutS = 60;
cfg.moon.url = 'https://ssd.jpl.nasa.gov/api/horizons.api';

%% 物理常量（通常无需修改）
cfg.constants.c = 299792458;
cfg.constants.omegaEarth = 7.2921151467e-5;
cfg.constants.muEarth = 3.986004415e14;
cfg.constants.muMoon = 4.902818954e12;
cfg.constants.earthRadiusM = 6378137;
cfg.constants.j2 = 1.08262668e-3;
end
