function cfg = lugreConfig()
% LUGRECONFIG mainEpoch（逐历元）与 mainBatch（多历元批处理）共用配置。
% 单位：距离 m，速度 m/s，钟差 c*dt (m)。
root = fileparts(mfilename('fullpath'));

%% 数据与运行范围
cfg.data.root = '/Users/ignat/Documents/LuGRE';
cfg.data.window = 'OP5_0';
cfg.data.sampleRateHz = 1;          % 用户确认整个数据集均为 1 Hz
cfg.data.rinex = 'auto';            % 'auto' 按观测 GPST 日期选 BRDC 文件；也可填具体路径
cfg.data.rinexDir = fullfile(root, 'data', 'rinex'); % auto 模式下 BRDC00IGS_R_YYYYDDD0000_01D_MN.rnx 所在目录
cfg.ephemeris.maxAgeS = 4*3600;     % 星历 Toe 与发射时刻最大间隔，超出的观测剔除
cfg.run.maxEpochs = Inf;            % 调试时可设 5、60；Inf 表示全部
cfg.run.skipEpochs = 0;             % 逐历元：跳过开头的 NAV 历元数（OP14 设 3，避开 178 s 空档）
cfg.run.splitGaps = false;          % 逐历元：NAV 有空档时按连续段分别处理（每段如独立运行，从 NAV 位置重新初始化）
cfg.run.plot = true;
cfg.run.progressEvery = 10;         % 逐历元：每多少历元打印进度，0 表示关闭
cfg.truth.enabled = true;
cfg.truth.file = fullfile(root, 'data', 'truth.xlsx');
cfg.truth.cacheFile = fullfile(root, 'results', 'cache', 'truth.mat'); % xlsx 解析结果缓存，源文件变化自动重建
cfg.truth.maxGapS = 1.5;            % 仅在连续真值样本间插值，不跨数据缺口

%% 逐历元方法开关（true 开启，false 关闭）
cfg.enable.dpe = true;
cfg.enable.ls = true;
cfg.enable.cn0Dpe = false;          % DPE 残差按 C/N0 加权
cfg.enable.cn0Ls = false;           % LS 残差平方按 C/N0 加权
cfg.enable.orbit = true;            % 轨道传播 + DPE 位置先验；LS 仍为纯观测对照
cfg.enable.j2 = false;              % 地球 J2 摄动；逐历元仅 orbit=true 时生效，批处理同样使用
cfg.enable.sagnac = true;           % 地球自转传播修正，所有方法共用
cfg.enable.legacyIonosphere = false;% 原经验公式不适合地月几何，仅供复现旧结果；批处理始终关闭

%% 逐历元初始状态（继续使用 NAV 初始位置、分段 NAV 速度）
cfg.initial.positionOffsetM = [0; 0; 0];
cfg.initial.clockBiasM = -140000;

%% 逐历元 DPE 搜索网格：相对上一有效解的偏移
cfg.dpe.switchEpoch = 201;
cfg.dpe.coarseOffsetsM = -20000:1000:20000;
cfg.dpe.fineOffsetsM = -10000:500:10000;
cfg.dpe.clockOffsetsM = -25000:1000:10000;
cfg.dpe.orbitWeight = 1e-5;         % 保留原值；距离惩罚系数，尚未标定
cfg.dpe.usePredictionCenter = true; % 有轨道预测时以预测位置为搜索中心
cfg.dpe.refineStepsM = [100, 10];   % 局部逐级精化；[] 关闭，不保证全局最优
cfg.dpe.rankTolerance = 1e-10;
cfg.dpe.chunkSize = 200000;         % 分块计算候选点，降低内存峰值

%% 逐历元 LS 迭代
cfg.ls.maxIterations = 30;
cfg.ls.positionToleranceM = 0.01;
cfg.ls.clockToleranceM = 0.01;
cfg.ls.rankTolerance = 1e-10;

%% 多历元批处理（mainBatch 使用；以 RAW 历元为轴，不用 NAV 初值和速度）
% dpeIrls: 加权 L1，IRLS 求解，样条钟差；ls: 同模型加权 L2 对照；
% dpeGrid: 6 维轨道网格 + 每历元中位数钟差，不用 IRLS。
cfg.batch.methods = {'dpeIrls','ls','dpeGrid'};
% 运动模型。'freeFlight'：tRef 的 ICRF 位置速度 + 地月引力积分（转移段）；
% 'lunarStatic'：月固系（DE405 主轴系，moonLibration）固定位置 p，r_ICRF = r_moon + R'*p，
% 用于着陆后的月面窗口（OP38 起）。该模型下 dpeGrid 为 3 维网格，先验盒只用 priorPositionM。
cfg.batch.dynamics = 'freeFlight';
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
cfg.batch.positionToleranceM = 1;     % 0.05 m 时 OP9 在 IRLS 容差量级来回抖动、报未收敛（解差 <2 m）
cfg.batch.irlsIterations = 100;       % L1 代价的迭代重加权次数上限
cfg.batch.irlsEpsilonM = 0.01;
cfg.batch.irlsToleranceM = 1e-4;
% dpeIrls 可选的最弱方向一维网格。消融（experiments/gridAblation.m，OP5/9/14 × 13 种初值偏差，
% 最大 1,000 km / 200 m/s）：开/关代价完全相同、解差 <2 m，关闭快 10~50 倍。线性化后 L1 代价是凸的。
cfg.batch.weakSearch.outerIterations = 0; % 前几次外迭代执行搜索；0 关闭，设 2 可作全局最优核验
cfg.batch.weakSearch.halfWidthM = 50000;  % 位置分量半宽
cfg.batch.weakSearch.stepM = 250;
% dpeGrid（solveBatchGrid）：白化坐标下逐级细化。白化后 1 个单位约相当于 1 m 伪距残差。
cfg.batch.grid.priorPositionM = 1e5;  % 首次外迭代的初值不确定盒：位置 ±100 km
cfg.batch.grid.priorVelocityMps = 100;%                       速度 ±100 m/s
cfg.batch.grid.pointsPerAxis = 5;     % 每轴点数（奇数），每级 5^6 = 15625 点
cfg.batch.grid.finalStep = 0.1;       % 白化坐标的最终步长
cfg.batch.grid.refineHalfWidth = 50;  % 之后外迭代（重新线性化）的白化半宽
cfg.batch.grid.maxLevels = 80;        % 每次外迭代的最大层数
cfg.batch.grid.chunkSize = 1024;      % 每批计算的网格点数，控制内存

%% 扩展卡尔曼滤波（mainFilter；移植自 EKF-LuGRE。状态 [r; v; 钟差; 钟漂]，伪距 + Doppler，RAW 历元为轴）
% ekfOrbit: 地心 ICRF 地球+月球动力学预测（含状态转移矩阵）；ekfCvBlend: ECEF 匀速 EKF 后与轨道预测按固定权重加权。
% 两种方法的噪声参数保持原项目的取值（原项目以 NAV 为参考调参，未按真值标定）。
cfg.filter.methods = {'ekfOrbit','ekfCvBlend'};
cfg.filter.minInitSatellites = 5;     % 初始化历元最少卫星数（伪距 LS 有 1 个多余观测）
cfg.filter.initIterations = 30;
cfg.filter.initStepToleranceM = 1e-4;
cfg.filter.initStd = [2000 2000 2000 20 20 20 2000 50]; % 初始标准差：位置 m、速度 m/s、钟差 m、钟漂 m/s
cfg.filter.initCovariance = 'ls';   % 'ls': 初始化 LS 的形式协方差；'fixed': 原项目的 initStd 对角阵（远距离过度自信，OP14/OP5 发散）
cfg.filter.gateSigma = 8;             % 归一化新息超过此值的观测剔除
cfg.filter.maxStepS = 5;              % 轨道 RK4 子步长上限
cfg.filter.cn0ReferenceDbHz = 35;     % 观测标准差 σ = σ_ref·10^((35-C/N0)/20)，再限幅
cfg.filter.rangeSigmaLimitsM = [5, 150];
cfg.filter.rangeRateSigmaLimitsMps = [0.08, 5];
cfg.filter.ekfOrbit.accelerationPsd = 1e-3;  % 未建模加速度 PSD (m/s^2)^2/Hz
cfg.filter.ekfOrbit.clockDriftPsd = 2;       % 钟漂随机游走 PSD (m/s^2)^2/Hz
cfg.filter.ekfOrbit.rangeSigmaM = 120;       % 35 dB-Hz 时伪距标准差
cfg.filter.ekfOrbit.rangeRateSigmaMps = 4;   % 35 dB-Hz 时变化率标准差
cfg.filter.ekfCvBlend.accelerationPsd = 9;
cfg.filter.ekfCvBlend.clockDriftPsd = 400;
cfg.filter.ekfCvBlend.rangeSigmaM = 8;
cfg.filter.ekfCvBlend.rangeRateSigmaMps = 0.15;
cfg.filter.ekfCvBlend.kalmanWeight = 0.7;    % 位置/速度 = w·EKF + (1-w)·轨道预测
% ekfLunarStatic：月面静止模型（月固系固定位置 p，无过程噪声；状态 [p; 钟差; 钟漂]），需在 methods 中显式加入。
% 钟与观测噪声沿用 ekfOrbit 的原项目取值，未标定。
cfg.filter.ekfLunarStatic.clockDriftPsd = 2;
cfg.filter.ekfLunarStatic.rangeSigmaM = 120;
cfg.filter.ekfLunarStatic.rangeRateSigmaMps = 4;

%% 观测筛选与信噪比权重
cfg.observation.minCount = 4;
cfg.observation.cn0ReferenceDbHz = 30;
cfg.observation.weightLimits = [0.01, 100]; % w=10^((C/N0-reference)/10)，再限幅
cfg.observation.legacyZenithDelayS = 5e-9;

%% 轨道传播（地球 + 月球 + 可选地球 J2）
cfg.orbit.resetIntervalS = 60;      % 逐历元：按实际时间重置，不按样本编号
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
