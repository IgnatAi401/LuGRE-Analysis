# LuGRE

## 本机运行

- 最终版入口：`main5.m`。
- 原始数据目录：`/Users/ignat/Documents/LuGRE`，遥测文件位于 `L0/TLM`。
- RINEX 文件保留在本项目根目录；入口脚本按自身位置加载依赖和星历。
- 当前入口固定分析 `OP1_0`，对应 `BRDC00IGS_R_20250150000_01D_MN.rnx`。切换观测窗口时也需要核对星历日期。
- `settings/config.txt` 的输出目录已改为项目下的 `results`；当前主脚本未使用此配置，后续写出文件时需要先创建目录。

运行：

```sh
/Applications/MATLAB_R2026b.app/bin/matlab -sd /Users/ignat/Developer/LuGRE-Analysis -batch "main5"
```

`main5` 的开关和参数统一在 [main5_config.m](main5_config.m)。轨道开启时需要 Aerospace Toolbox，并联网请求 JPL Horizons；关闭轨道后不请求月球星历。真值对比也需要 Aerospace Toolbox 进行坐标转换。

### main5 开关

直接编辑 `main5_config.m` 中的 `true` / `false`：

| 配置 | 作用 | 默认 |
| --- | --- | --- |
| `cfg.enable.dpe` | 网格 DPE 定位 | 开 |
| `cfg.enable.ls` | 最小二乘定位，对照方法独立迭代 | 开 |
| `cfg.enable.cn0Dpe` | DPE 的 C/N0 权重 | 关 |
| `cfg.enable.cn0Ls` | LS 的 C/N0 权重 | 关 |
| `cfg.enable.orbit` | 轨道传播及 DPE 位置先验 | 开 |
| `cfg.enable.j2` | 轨道中的地球 J2 摄动 | 开 |
| `cfg.enable.sagnac` | 传播期间地球自转修正 | 开 |
| `cfg.enable.legacyIonosphere` | 旧经验电离层项，仅用于复现 | 关 |

DPE 与 LS 至少开启一个。两者同时开启时，轨道位置反馈使用 DPE；只开 LS 时使用 LS。LS 始终是纯观测解，不叠加轨道先验；`orbit=true` 时即使只开 LS，也会输出单独的轨道预测。两种方法仍保留 NAV 初始位置和轨道重置时的 NAV 速度辅助。

信噪比开关使用 C/N0（载噪密度比，单位 dB-Hz），权重为 `10^((CN0-reference)/10)`，按配置限幅；DPE 加权绝对残差，LS 加权残差平方。DPE 固定位置后用加权中位数找到最优离散钟差的相邻候选，避免重复计算每个钟差点的几何距离。默认关闭旧电离层经验项，因为它的角度定义不适用于此处几何；需要复现旧观测模型时可显式开启。

### 常用参数与输出

- `cfg.run.maxEpochs`：先设为 `5` 或 `65` 检查运行，`Inf` 处理所有历元。
- `cfg.dpe.*`：位置/钟差网格、切换历元、轨道先验权重、计算分块大小。`usePredictionCenter` 默认用轨道预测作为中心；`refineStepsM=[100,10]` 局部精化位置与钟差，设为 `[]` 关闭。
- `cfg.ls.*`：迭代次数、位置和钟差收敛阈值。
- `cfg.orbit.*`：按秒计的重置间隔与积分容差。
- `cfg.frame.*`：GPS-UTC 闰秒差、地球定向参数来源。默认使用 MATLAB IERS 数据；手动模式的零 UT1/极移值仅是近似，不会自动作为失败回退。

`results.dpe`、`results.ls` 包含位置、以米表示的钟差、有效性标记和与 NAV 的距离差；`results.orbit` 为轨道预测。没有开启的方法保持 NaN，不绘制假对照。DPE 还输出 `onBoundary`，用于识别网格边界解。`main5` 保留 `pos/pos2/dt/dt2/dd/dd2/dd3/NAV` 等常用工作区变量。

也可以从 MATLAB 调用函数，自定义一次运行而不修改配置文件：

```matlab
cfg = main5_config();
cfg.run.maxEpochs = 65;
cfg.run.plot = false;
cfg.enable.orbit = false;
cfg.enable.cn0Dpe = true;
results = runLugreAnalysis(cfg);
```

整个数据集采样率由用户确认为 **1 Hz**，入口核查所选 NAV 序列的秒间隔。每颗卫星只保留一个信号（最小 signalId 优先，同信号取最高 C/N0），以星座和 PRN 去重；Toc 用于卫星钟差，Toe 用于轨道传播并处理跨周。信号群延迟与系统间偏差尚未完整建模。

### 真值对比

`main5_config.m` 的 `cfg.truth.enabled` 默认开启。`temp.xlsx` 的 C 列为 GPS 秒，F:H 为地心 J2000 位置（km）；按 GPS 秒线性插值，再使用同一套 IAU/IERS 转换得到 ITRF 米。禁止外推或跨越超过 `cfg.truth.maxGapS` 的缺口。真值只用于事后评估，不参与算法。`results.truth` 保存真值位置、有效性及 DPE、LS、轨道预测、NAV 的 XYZ 误差、三维距离误差、RMSE 和中位数。真值转换需要 Aerospace Toolbox，即使关闭轨道也一样；不需要对比时可关闭 `cfg.truth.enabled`。

历史结果重画：`replotLugreTruth('try1.mat')`，无需重跑定位；输出到 `results/try1_truth/`，包含三张 `.fig`、对应 PNG 和带真值统计的 `truth_comparison.mat`，原 MAT 不变。

### 坐标与时间约定

月球星历使用地心 ICRF 赤道轴、几何矢量和 UTC 时间（Horizons `TIME_TYPE=UT`）。GPST 先减去配置中的 GPS-UTC；星历按实际观测时间插值。只把位置和速度从 km、km/s 转为 m、m/s，JD 时间列不变。航天器使用 IAU 2000/2006 的 ICRF/GCRF 与 ITRF 转换，含岁差、章动、地球自转和 IERS 极移；速度使用旋转矩阵导数，J2 在地固系计算后转回惯性系。

新 `main5` 使用 `lugreDpe/lugreLs/lugreDynamics`；旧版 `DPE_rinex_Moon_orbit3`、`leastSquarePos` 等保留供历史入口使用。旧入口不读取这份新配置中的方法开关。

离线回归检查（坐标测试需要 Aerospace Toolbox）：

```matlab
addpath('tests');
test_main5
```

在线集成检查：`test_main5_integration`。需要本机 LuGRE 数据、Aerospace Toolbox 和网络；会检查默认网格的 65 个历元、轨道重置、六种求解/加权组合及绘图。

本轮验证：离线回归（含与四维穷举等价性）通过；轨道开启的前 65 个历元 DPE/LS 均为 65/65 有效，六种求解/加权组合、J2 开关和绘图检查通过。关闭轨道的完整 OP1_0 共 1,387 个历元，两种方法均为 1,387/1,387 有效，末历元已计算；与 NAV 距离差中位数分别为 DPE 2,582.622 m、LS 386.834 m。DPE 有 1,170 个历元在至少一级搜索中触及边界，不能据此声称充分收敛或精度提升。全序列开启轨道的检查未完成，轨道权重尚未标定。

初次算法审查及其局限见 [MAIN5_REVIEW.md](MAIN5_REVIEW.md)。该文件中的原始代码行号仅对应重构前版本。

This repository contains **MATLAB analysis scripts** for the [LuGRE](https://etd.gsfc.nasa.gov/our-work/lunar-gnss-receiver-experiment-lugre/) (Lunar GNSS Receiver Experiment) [Mission Data](https://zenodo.org/records/16411687).

## Overview

- **Mission:** LuGRE — Lunar GNSS Receiver Experiment
- **Data:** Public LuGRE/BGM1 mission data, including raw GPS/Galileo L1/E1 and L5/E5 observables and signal sample batches
- **Language:** MATLAB
- **Purpose:** GNSS-based positioning, navigation, and orbit determination in cislunar and lunar scenarios

## Data Sources and Ephemerides

The analysis uses:

- **NASA hybrid broadcast ephemeris** for GNSS satellite orbits
- **JPL Horizons lunar ephemeris** for lunar position and geometry

These inputs support the reconstruction and analysis of LuGRE observations in the Moon-transfer and lunar environment.

## Algorithms and Methods

The implemented workflow includes:

- **Observation-level Direct Position Estimation (DPE)**  
  Direct estimation of position from GNSS observables at the observation level.

- **Least-squares estimation**  
  Used for parameter estimation and navigation solution computation.

- **Restricted three-body orbital dynamics optimization**  
  Used for improved orbit determination and dynamical consistency in the Earth–Moon system.

## Key Highlights

- Integrates **NASA hybrid broadcast ephemeris** and **JPL Horizons lunar ephemeris**
- Combines **observation-level DPE**, **least-squares estimation**, and **restricted three-body orbital dynamics optimization**
- Supports **GNSS-based cislunar/lunar navigation and orbit-determination analysis**
- Built on **public LuGRE mission data** with **MATLAB** analysis scripts
