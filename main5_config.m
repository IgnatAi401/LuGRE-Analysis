function cfg = main5_config()
% MAIN5_CONFIG main5 的统一配置。距离 m，速度 m/s，钟差 c*dt (m)。
root = fileparts(mfilename('fullpath'));

%% 数据与运行范围
cfg.data.root = '/Users/ignat/Documents/LuGRE';
cfg.data.window = 'OP5_0';
cfg.data.sampleRateHz = 1;          % 用户确认整个数据集均为 1 Hz
cfg.data.rinex = 'auto';            % 'auto' 按观测 GPST 日期选 BRDC 文件；也可填具体路径
cfg.data.rinexDir = root;           % auto 模式下 BRDC00IGS_R_YYYYDDD0000_01D_MN.rnx 所在目录
cfg.ephemeris.maxAgeS = 4*3600;     % 星历 Toe 与发射时刻最大间隔，超出的观测剔除
cfg.run.maxEpochs = Inf;             % 调试时可设 5、60；Inf 表示全部
cfg.run.plot = true;
cfg.run.progressEvery = 10;         % 每多少历元打印进度，0 表示关闭
cfg.truth.enabled = true;
cfg.truth.file = fullfile(root, 'temp.xlsx');
cfg.truth.maxGapS = 1.5;            % 仅在连续真值样本间插值，不跨数据缺口

%% 方法开关（true 开启，false 关闭）
cfg.enable.dpe = true;
cfg.enable.ls = true;
cfg.enable.cn0Dpe = false;           % DPE 残差按 C/N0 加权
cfg.enable.cn0Ls = false;            % LS 残差平方按 C/N0 加权
cfg.enable.orbit = true;            % 轨道传播 + DPE 位置先验；LS 仍为纯观测对照
cfg.enable.j2 = false;               % 仅 orbit=true 时生效
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
cfg.dpe.usePredictionCenter = true; % 有轨道预测时以预测位置为搜索中心
cfg.dpe.refineStepsM = [100, 10];   % 局部逐级精化；[] 关闭，不保证全局最优
cfg.dpe.rankTolerance = 1e-10;
cfg.dpe.chunkSize = 200000;         % 分块计算候选点，降低内存峰值

%% LS 迭代
cfg.ls.maxIterations = 30;
cfg.ls.positionToleranceM = 0.01;
cfg.ls.clockToleranceM = 0.01;
cfg.ls.rankTolerance = 1e-10;

%% 多历元批处理 DPE（mainBatch 使用；以 RAW 历元为轴，不用 NAV 初值和速度）
cfg.batch.methods = {'dpe','ls'};     % dpe: 加权 L1 + 最弱方向网格；ls: 同模型加权 L2 对照
cfg.batch.clockModel = 'spline';      % 'spline': 三次 B 样条；'polynomial': 整段多项式
cfg.batch.clockKnotSpacingS = 120;    % spline 节点间距；OP5 实测 120 s 内钟差近似二次
cfg.batch.clockDegree = 2;            % polynomial 模式的阶数
% 信号偏差：OP5 每种信号只来自 1~2 颗星，偏差参数等同于单星偏差，会吸收区分径向与钟差的
% 视线差异信息（实测开启后径向误差约 1,000 km、关闭约 10~30 km）。默认关闭。
cfg.batch.signalBias = false;         % 估计各 signalId 相对最常用信号的常值偏差
cfg.batch.minSignalObservations = 20; % 少于此数的信号不参与，避免偏差不可估
cfg.batch.cn0Weight = false;          % 按 C/N0 加权，沿用 observation 的参考值与限幅
cfg.batch.initMinSatellites = 4;      % 快照 LS 初值所需最少卫星数
cfg.batch.initMinEpochs = 10;         % 快照解少于此数时无法初始化
cfg.batch.allowNavInit = false;       % 快照不足时用 NAV 位置初始化（属 NAV 辅助，需注明）
cfg.batch.maxOuterIterations = 10;    % 轨道重积分 + 重新线性化
cfg.batch.positionToleranceM = 0.05;
cfg.batch.irlsIterations = 100;       % L1 代价的迭代重加权次数上限
cfg.batch.irlsEpsilonM = 0.01;
cfg.batch.irlsToleranceM = 1e-4;
cfg.batch.gridHalfWidthM = 50000;     % 最弱方向一维网格：位置分量半宽
cfg.batch.gridStepM = 250;
cfg.batch.gridOuterIterations = 2;    % 前几次外迭代执行网格搜索

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
