# 比较1正式 Monte Carlo

本目录代码只经过文本级静态检查，尚未由 MATLAB 解析或执行，也未运行验证或仿真。代码准备阶段没有启动 MATLAB。

固定实验：SNR `-20:5:30` dB，每个 SNR 300 次，种子 `20261001:20261300`。相同 trial 在各 SNR 下重新设置相同 seed，以复用 RIS 码字与标准化噪声；trial 间不传播任何状态。

## 正式运行及断点续跑

首次运行和断点续跑使用相同命令：

```matlab
root = 'E:\file\matlab\MultiRIS3D_paper_baseline_A_20260918_173645';
addpath(fullfile(root,'formal'));
run_formal_monte_carlo
```

入口会恢复调用前的 MATLAB path 与 RNG，不调用 `savepath`。每个完成的 trial 都原子更新对应 SNR 检查点。若配置、源码哈希、SNR、trial 数或 seed 不匹配，续跑会明确报错。

检查点保存于：

```text
formal_results/checkpoints/checkpoint_snr_m020.mat
...
formal_results/checkpoints/checkpoint_snr_p030.mat
```

## 只汇总

```matlab
root = 'E:\file\matlab\MultiRIS3D_paper_baseline_A_20260918_173645';
addpath(fullfile(root,'formal'));
summary = summarize_formal_monte_carlo();
```

汇总仅读取现有检查点，不调用信号生成、信道估计或位置推断。允许对未完成检查点生成阶段性汇总；每项统计都记录实际样本数。输出：

```text
formal_results/summary/formal_summary.mat
formal_results/summary/formal_summary.csv
```

## 只绘图

```matlab
root = 'E:\file\matlab\MultiRIS3D_paper_baseline_A_20260918_173645';
addpath(fullfile(root,'formal'));
summaryFile = fullfile(root,'formal_results','summary','formal_summary.mat');
figureDir = fullfile(root,'formal_results','summary','figures');
outputFiles = plot_formal_monte_carlo(summaryFile, figureDir);
```

绘图只加载现有 `formal_summary.mat`。输出 RMSE、失败率、运行时间和 Joint 误差分布四组 PNG/PDF/FIG。

## 统计口径

- 粗定位指标使用具有有限粗定位输出的样本。
- 精修指标使用具有有限精修输出的样本，包括有限但未正常收敛的输出。
- 无输出不填零。
- 无输出率和有限未收敛率以已完成 trial 为分母；另保存“有限未收敛/有限精修输出”的条件比例。
- p95 使用排序后的 nearest-rank 定义 `ceil(0.95*N)`。
- 正式方法输出固定为精修输出，不根据真实误差在粗定位和精修间择优。
