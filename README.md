# 小论文比较1：Multi-RIS 3D 方案 A

**Multi-RIS 3D adapted to the single-frame two-user scenario, with oracle support and a truth-centered coarse search region.**

本目录是独立、单帧、双静止用户的方案 A 实验。它不调用 VB，不生成轨迹，不使用位置先验惩罚，不实现比较2，也不运行正式 Monte Carlo。

## 来源与授权改动

- 原仓库：`E:\file\matlab\Tracking_code`
- 历史实验快照：`327664ee229c7680f9a8bd71b008ff14ba2a0b48`
- `src/` 下的十个文件由该提交直接抽取，内容未修改。
- `simulation_config.m` 保留历史原文；`paper_baseline_A_config.m` 核对固定物理量后，将历史 `TFT` 明确覆盖为实验要求的 `TTT`。
- 新增的 `separate_noiseless_observation.m` 调用原始物理生成器，并逐行复用历史 DFT 分离与归一化步骤。
- 方案 A runner 将原距离域位置精修替换为已有的参数域单位权重 `refine_position_mle_eq57_identity.m`；其残差顺序为 `[xi; zeta; tau_seconds]`。

完整来源、blob、首次引入提交和 SHA-256 见 `provenance/source_manifest.csv`；修改说明见 `provenance/modification_manifest.md`。

## MATLAB 命令

静态检查：

```matlab
root = 'E:\file\matlab\MultiRIS3D_paper_baseline_A_20260918_173645';
oldPath = path;
oldRng = rng;
cleanup = onCleanup(@() restoreSession(oldPath, oldRng));
addpath(fullfile(root,'experiment'));
report = static_check()

function restoreSession(oldPath, oldRng)
path(oldPath);
rng(oldRng);
end
```

运行限定的 A/B/C 验证：

```matlab
root = 'E:\file\matlab\MultiRIS3D_paper_baseline_A_20260918_173645';
oldPath = path;
oldRng = rng;
cleanup = onCleanup(@() restoreSession(oldPath, oldRng));
addpath(fullfile(root,'validation'));
validation = run_all_validations();

function restoreSession(oldPath, oldRng)
path(oldPath);
rng(oldRng);
end
```

入口自身也使用 `onCleanup` 恢复 MATLAB path 和 RNG；不调用 `savepath`。

加载最终结果：

```matlab
root = 'E:\file\matlab\MultiRIS3D_paper_baseline_A_20260918_173645';
saved = load(fullfile(root,'results','validation_all.mat'),'validation');
validation = saved.validation;
summary = struct2table(validation.C.summary)
```

单样本结果位于 `results/trial_seed*_snr*dB.mat`，C 层汇总位于 `results/validation_C_summary.csv`。

## 算法链

`generate_noisy_separated_signal` → `estimate_delay_coarse` → `estimate_spatial_coarse` → `refine_channel_mle` → `estimate_position_coarse` → `refine_position_mle_eq57_identity`。

粗搜索范围为 UE1 真值每轴 ±1 m，步长 0.2 m，共 `11^3=1331` 点。真值只在实验驱动层构造搜索边界，不参与候选选解或连续精修。

## 已执行验证

- A：真实传播参数几何一致性和独立有限差分 Jacobian 检查，通过。
- B：共享生成器的严格无噪声完整流程，通过；没有以高 SNR 代替无噪声。
- C：SNR `[0,30]` dB、种子 `20260918:20260922`，恰好 10 次；全部保存且没有剔除。

数值摘要见 `results/validation_report.txt`，完整结构见 `results/validation_all.mat`。
