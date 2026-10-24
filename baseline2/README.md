# Baseline 2：Prior-aided Multi-RIS 3D

本目录独立实现 Baseline 2 的信道参数 Fisher 信息模块。Baseline 1 的源码、配置、验证及已有 3300 次实验结果均不修改，也不复制 `src/`。本轮不实现完整 Baseline 2。

## 研究目标与计划流程

保留 Baseline 1 的时延/空间频率粗估计和信道 MLE，后续引入与 Structured VB 相同的六维高斯位置先验，实施先验辅助位置初始化与 Fisher 加权位置 MAP。区别在第二阶段的位置推断；现有 Baseline 1 不会被替换或重跑。

六维位置顺序为 `[pT; pR]`。每个 trial 的先验均值、协方差必须与 Structured VB **完全相同**，由统一实验入口提供，不得由 Baseline 2 独立生成。按本轮任务提供的信息，当前 Structured VB 研究分支使用 `C0 = blkdiag(0.5*eye(3), 100*eye(3))`；本模块不硬编码或使用该协方差。

## 已实现：单路径等效 Fisher 信息

入口：`inference/compute_channel_fim.m`。

```matlab
result = compute_channel_fim(qHat, rhoHat, W, separatedNoiseVariance, cfg);
```

`qHat` 为实数 `3×1` 列向量 `[xi; zeta; nu]`，其中 `nu = cfg.ofdm.subcarrier_spacing * tau`，不是以秒为单位的时延。`W` 为本帧该 RIS 的 `N×G` 基础反射矩阵。噪声模型为 `CN(0, sigmaSep2*I)`，输入方差为 DFT 分离后每个复样本的 `E|n|²`，不是单个实部的方差，也不是分离前方差；现有生成器分离后方差等于分离前方差除以 DFT 块数。本函数仅接受正噪声方差。

信号原子与现有 MLE 一致：

\[
\mu=(W^{\mathsf T}a)\otimes d,\qquad
a=a_x\otimes a_y,\qquad d_k=e^{-j2\pi k\nu},\quad k=0,\ldots,K-1.
\]

阵列索引从零开始，x 轴慢变、y 轴快变；使用 `W.'` 非共轭转置。解析导数使用 `a_x' = j(2*pi/lambda)*dx*n_x.*a_x`、对应的 y 轴导数及 `d_nu = -j*2*pi*k.*d`，再通过同一 `W.'` 与 Kronecker 顺序构造 `D=[mu_xi,mu_zeta,mu_nu]`。这里撇号表示求导，不是 MATLAB 转置。

消去未知复增益的等效信息为：

\[
F=\frac{2|\widehat\rho|^2}{\sigma_{\rm sep}^2}
\operatorname{Re}\left[D^{\mathsf H}D-
\frac{(D^{\mathsf H}\mu)(\mu^{\mathsf H}D)}{\mu^{\mathsf H}\mu}\right].
\]

这相当于对未知复增益的两个实参数作 Schur 消元；不能省略第二项。实现不构造大型投影矩阵，仅对最后的实 `3×3` 矩阵做对称化，不求逆、不截断特征值、不正则化。返回 `fim`、`eigenvalues`、`conditionNumber`、`parameterOrder`，以及 `psdTolerance`、`isPositiveSemidefinite`、`isNearlySingular` 诊断。PSD 容差为相对浮点舍入容差；近奇异判据为 `rcond(F)<sqrt(eps)`，只记录、不调整参数。可选第二输出 `[result, model]` 提供原子与解析导数，供轻量验证使用。

## 与现有估计器的衔接

后续由统一入口调用现有 `src/baselines/multiris3d/estimate_channel_parameters.m`，本模块不重复实现估计器。对其结果 `channel` 的第 `j` 条活动路径：

```matlab
ell = channel.activeRisIndices(j);
qHat = [channel.xiEstimate(j); channel.zetaEstimate(j); ...
    cfg.ofdm.subcarrier_spacing * channel.delayEstimate(j)];
info = compute_channel_fim(qHat, channel.rhoEstimate(j), ...
    baseReflections(:,:,ell), channel.separatedNoiseVariance, cfg);
```

Fisher 信息是估计点处的局部 plug-in 信息，不是低 SNR 下估计误差分布已被验证为高斯的保证。后续位置残差必须使用同样的 `[xi; zeta; nu]` 坐标，并处理周期分支。若改用 `r=[xi;zeta;tau]`，令 `A=diag([1,1,df])`，则 `F_r=A.'*F_q*A`，不可只改残差单位而保持权重不变。

## 最小验证与当前状态

运行 `validation/validate_channel_fim.m`：只读取已有 `results/trial_seed20260918_snr+30dB.mat` 中的配置、基础反射矩阵及信道估计结果，不加载观测数组，不生成信号、不调用估计器、不进行位置优化、不写文件。使用现有 `upa_steering.m` 与 MLE 的 `M(:)` 构造作独立中心差分参考，检查三个实际 RIS 估计点的解析导数、Fisher 维度、实对称性、半正定性及参数顺序，并报告特征值和条件数。

当前状态：2026-10-24 已运行一次轻量 MATLAB 验证，三个 RIS 的解析导数、维度、实对称性、半正定性和参数顺序检查全部通过。读取样本的分离后复噪声方差为 `3.73734753482e-19`。具体结果如下（特征值按升序）：

| RIS | ξ 导数相对误差 | ζ 导数相对误差 | ν 导数相对误差 | Fisher 特征值 | 条件数 |
| --- | --- | --- | --- | --- | --- |
| 1 | 1.390e-10 | 7.872e-10 | 1.416e-10 | 1.454714e8, 1.800447e8, 1.021654e10 | 70.2306 |
| 2 | 1.314e-10 | 7.725e-10 | 1.476e-10 | 9.143742e7, 1.356291e8, 1.156410e10 | 126.470 |
| 3 | 6.654e-10 | 1.924e-10 | 1.190e-10 | 3.273409e8, 5.169961e8, 2.873354e10 | 87.7786 |

本样本均为正定、非近奇异，没有进行正则化；条件数仅说明当前参数坐标下的本样本情况，不代表所有 trial。未进行定位或 Monte Carlo 验证。

## 后续待实现（本轮不实施）

1. 从统一入口接收同 trial 的六维位置先验，实施先验辅助初始化。
2. 在一致的归一化参数及周期分支下构建 Fisher 加权六维位置 MAP，结合先验信息；明确退化信息的处理策略。

目前只有 Fisher 模块，不能将本目录视作已完成或已验证的定位 Baseline。
