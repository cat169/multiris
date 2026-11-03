# Baseline 1-PAI 源码审查材料

本分支只发布已验证源码和本文档，供正式集成前审查；未修改算法、原 B1 或 Tracking 实验入口，未重新运行实验，也不合并 main。B1-PAI 是 Multi-RIS 3D with Projection-Assisted Initialization 数值增强实现，不是原论文未经修改的算法或 AML 直接复现。采用为论文 Baseline 的决定不替代下述数值风险披露。

## 1. 接口与集成适配点

```matlab
result = infer_baseline1_frame_projection_assisted( ...
    YSeparated, separatedNoiseVariance, baseReflections, cfg);
% 配对验证可选：复用对同一观测、同一配置计算的原 MLE。
result = infer_baseline1_frame_projection_assisted( ...
    YSeparated, separatedNoiseVariance, baseReflections, cfg, channel);
```

- `YSeparated`：子载波 × 基础反射配置 × RIS 的复分离观测。
- `separatedNoiseVariance`：原分离噪声方差；原信道接口保留该标量。
- `baseReflections`：RIS 元素 × 基础反射配置 × RIS。
- `cfg`：原配置结构；包含 `ris`、`ofdm`、`constants`、RIS 中心/局部轴/法向，以及 `paperBaselineA.oracleSupport/lowerBound/upperBound/gridStep`。当前实现仍要求 `geometry.static_user/mobile_user` 字段存在，以便移除后向初始化计算传递无用户真值的几何结构；不使用其值选候选。
- 可选 `channel`：原 `estimate_channel_parameters` 输出。接口检查标量结构、必要字段、支持集合及其顺序、向量长度、有限实数及可选噪声方差匹配；这些结构检查不证明观测同源性，仍由调用方保证，不重新执行MLE。

与原 `experiment/run_paper_baseline_A.m` 相同的主要字段为 `staticPositionEstimate`、`mobilePositionEstimate`、`activeRisIndices`、`searchInfo`、`channel`、`coarse`、`refined`、`label`。新增：

| 字段 | 用途 |
| --- | --- |
| `finite`、`exceptionIdentifier/Message` | 结果是否有数值有限且原目标可计算的候选；无可用精修时返回 NaN 位置及 `baseline1PAI:NoValidRefinement` |
| `fallbackTriggered`、`initializationSource` | 原几何路径或投影辅助路径 |
| `originalExceptionIdentifier/Message` | 保留原 NoValidCandidate 证据 |
| `candidates`、`selectedIndex`、`candidateCount` | 全部实际精修候选、原目标最小候选索引；正常路径有1个候选 |
| `initializationDiagnostics` | 投影、切向法向、交会秩/条件数、负射线参数、直线最近距离等逐网格诊断 |
| `initialObjective`、`finalObjective` | 原秒单位、单位权重完整位置目标 |
| `qualityFailed`、`qualityReasons`、`strictGradient` | 质量报告；不筛除有限位置、不参与选优 |
| `timing`、`extraPositionSolves`、`channelReused` | 阶段计时、额外精修次数、是否复用外部 MLE |

未来 Tracking 适配需选择新接口并更新方法标签，保留有限性、质量及回退统计。应先检查 `finite`/异常和 `refined` 是否存在，再读取 `refined.exitflag`：全部候选无效时 `refined` 是空结构；回退时 `coarse` 也是空结构，不能假设原粗定位输出存在。原异常仍应由实验层捕获。不得重复调用旧 B1 再让新接口重估相同信道；配对时可共用原 MLE/正常计算路径。本分支没有实施任何 Tracking 适配。

使用时只需将仓库 `src/baselines/multiris3d` 与 `src/structured/current/signal` 加入 MATLAB path；`private` 由 MATLAB 的父目录可见规则解析，不需要独立加路径。需要 Optimization Toolbox 的 `lsqnonlin`；其余优化使用原 `fminsearch`。无需 Tracking 的 analysis、reports、临时目录或数据文件作为推断运行时依赖。

## 2. 实际执行顺序

1. 未传第5参数或其为空时，调用原 `estimate_channel_parameters` 一次：1024点 IFFT 时延粗搜 → 二维空间粗搜 → 单初值 `fminsearch` 联合信道精化。传入有效 `channel` 时跳过该步。
2. 从原 `cfg.paperBaselineA` 提取支持和 Tx 搜索条件，调用原 `estimate_position_coarse` 一次。只捕获完整标识 `estimate_position_coarse:NoValidCandidate`，其他异常原样抛出。
3. **正常路径**：原几何初值 → 原 `refine_position_mle_eq57_identity` 一次 → 独立复核原目标/最近RIS距离 → 附加质量和时序。原估计位置及 `refined` 不变，无投影或额外精修。
4. **回退路径**：`b1_pai_grid` 根据原边界/步长生成网格，当前验证配置为1331点，未将1331写死。逐 Tx 网格、逐 active RIS 减去 Tx 方向的面内分量，对越界二维向量除以其模长；可行向量原值保留。正法向分量为 `sqrt(max(0,1-sum(ab.^2)))`，无内缩或法向裕量。
5. 使用配置局部轴/法向、原 unordered RIS pair、`[uI,-uJ] \\ (rj-ri)`、line-i 最近点、多 pair 位置等权融合。Tx网格与任一有效RIS中心完全重合时跳过并标记`txAtRisCenter`，不增加非零距离阈值。秩只计算一次并复用原`<2`判定。秩不足、非有限或原残差不可计算时不纳入初值；负射线参数和高条件数仅诊断，不新增约束。
6. 每个可计算初值评价原完整目标，按目标、网格索引、RIS索引升序。六维欧氏距离 `<=1e-6 m` 去重，最多保留8个初值。没有使用真实误差、历史 Trial 编号、测试种子、Proposed 或 B2 输出。
7. 每个保留初值调用同一个原位置精修函数。目标仍是每RIS的 ξ、ζ 与以秒计的 τ 残差平方和，单位权重、无先验/Fisher正则；原 trust-region-reflective 参数不变（500次迭代、5000函数评价、Function/StepTolerance=1e-12、OptimalityTolerance=1e-10）。
8. 独立复核初/终原目标。在一致性断言前检查六维位置、优化器及独立初/终目标的有限实数性；非有限候选原值保留、标为无效并继续。有限目标真实不一致仍抛原`fallback:*ObjectiveMismatch`。仅对3个明确的候选数值异常标识记录后继续，其他错误原样抛出；详见修复报告。原目标最小者胜出，精确平局保留较早候选。exitflag/一阶指标/质量标签都不用于选择。
9. 逐原因标记非正或非有限exitflag、一阶指标非有限/>1e-6、距有效RIS中心<=1e-6m、目标增加>1e-10*max(1,J0)、非有限位置/目标；严格一阶标志为<=1e-10。`firstOrderOptimality` 是原求解器 `output.firstorderopt`，不是保存的完整梯度向量。

## 3. 依赖和开销审查

新增并提交的5个 private helper：

- `b1_pai_build_initializations.m`：投影、RIS交会、原目标评价及诊断。
- `b1_pai_grid.m`：原 Tx 网格。
- `b1_pai_rank_initializations.m`：确定排序、去重、前8个。
- `b1_pai_evaluate_candidates.m`：原精修、目标复核、有限候选选优。
- `b1_pai_objective.m`：原残差表达式及分项/最近RIS距离。

已在基准提交中存在、未修改且不重复提交的必要依赖：`estimate_channel_parameters.m`、`estimate_delay_coarse.m`、`estimate_spatial_coarse.m`、`refine_channel_mle.m`、`estimate_position_coarse.m`、`refine_position_mle_eq57_identity.m`，均位于 `src/baselines/multiris3d`；间接依赖 `src/structured/current/signal/upa_steering.m`。验证时解析的 Tracking steering 与本仓库 steering 文本一致（仅文件换行字节不同），无需复制 Tracking 文件。

| 检查点 | 源码发现（未优化） |
| --- | --- |
| 重复 MLE | 单次接口默认1次；第5参数复用时0次。外层若先跑旧接口再调用四参数新接口，会造成调用方重复，集成时需避免 |
| 重复原粗定位 | 每次接口1次；回退是在失败后另做投影网格扫描，不重复调用原粗定位 |
| 目标求值 | 正常路径在原精修后复核终点目标；回退对全部可计算网格初值求原目标，8个精修内部再次评价，返回后再复核。部分计算重复，但复核服务于排序/一致性诊断，未删除 |
| 相同初值 | 前8个按1e-6m六维距离去重，不重复完全相同初值；未强制起点大范围分散，不同初值仍可收敛到同一终点 |
| 矩阵开销 | 主要是原信道二维搜索、1331网格及最多8次精修；pair矩阵仅3×2。重复`rank(A)`已消除；`cond(A'*A)`和逐网格诊断保留，未发现新增的大型稠密矩阵求逆 |
| 磁盘 | 推断接口和 private helper 无 load/save/fopen 或其他磁盘读写；checkpoint 属于外层验证程序，未随本分支发布 |

## 4. 已有验证与对应配置

依据既有 `BASELINE1_PAI_IMPLEMENTATION_VALIDATION.md`、`historical_regression.csv`、`high_snr_regression.csv`、`support_snr_summary.csv` 和源码指纹清单；本轮只读引用，没有重跑或把报告作为代码运行时输入。验证数据未上传本分支。

- MATLAB R2022a；原G0布局、40子载波、2.5MHz间隔；原噪声/111参考幅度不变，双路径不额外功率归一化。4I实验条件不向B1目标增加高斯先验。
- 历史162个失败全部复现：选中候选编号、exitflag、迭代数一致，位置最大差0m；Joint RMSE=61.357790476019m。历史40个正常Trial通过，无回退、位置逐位相同。
- 60个中高SNR正常回归（4状态×0/10/20dB×5）全部通过，无回退、原粗定位一致、最终位置一致。
- 独立800例为4状态×−15/−10dB×100；每Trial原B1与PAI共享观测、反射矩阵和一次原MLE。种子不参与推断源码；种子映射/冲突核对及完整数据保留在既有本地验证材料中。
- 原B1有限777/800，PAI有限800/800；新增23例全部>10m，其中18例质量不合格、2例近RIS中心、3例>100m。1m/10m定位成功率两者均为44.5%/92.5%（分母800）。所以不能以有限率100%声称精度提升或普遍数值可靠。
- 本机固定工作负载平均原B1约0.352s、PAI约0.382s；每Trial平均额外精修0.23次。实际触发回退时为8次，不应将总体均摊开销当作单次回退开销。
- 两者继承原B1的 oracle Tx ±1m搜索区域；无额外Rx真值辅助或误差选优。有限极端错解/未收敛结果保留，质量标签只报告。

## 5. 进一步审查项（本轮不修复）

1. 返回契约：无可用回退候选返回NaN及空`refined`；未来调用方必须适配，不能沿用无条件访问`refined.exitflag`的写法。其他原异常仍会抛出，需要实验层独立捕获。
2. 非有限精修结果导致断言中断的问题已修复并经隔离替身验证；有限目标不一致、维度或未知程序错误仍故意立即抛出。数值异常白名单保守，不按错误消息猜测未知错误可恢复。
3. 可选外部MLE的明显结构错误已检查，但不能证明其与观测严格同源。`geometry.static_user/mobile_user` 字段虽不用于选择，现有配置契约仍要求字段存在。
4. 圆盘边界产生切向方向；负射线和近RIS中心结果未被过滤。这与已验证规则一致，任何修改都须另行验证，不能在本次发布中悄悄调整。
5. 原秒单位无权重时延项很弱，原目标最低和正exitflag都不保证真实定位准确。原始/PAI的有限分母、质量和固定误差门限成功率必须分别报告。

## 6. 源码身份

原B1基准提交：`431338bf92312c89eb76787adfffb44ceffbcd74`；初版B1-PAI验证/审查提交：`fe8253e8f9885af896cf7c037b2fc86427194c14`。本次最小工程修复只涉及入口和2个private helper，262例逐候选回归及23项边界测试通过，详见[工程修复报告](BASELINE1_PAI_HARDENING_REPORT.md)。下表为当前工作文件原始字节SHA-256，修复前后对照见报告。既有依赖未修改，换行字节可能受Git checkout设置影响。

| 文件（相对仓库根） | 验证 SHA-256 |
| --- | --- |
| `src/baselines/multiris3d/infer_baseline1_frame_projection_assisted.m` | `06fe32fc43d79c3f4ca4192e529b1a8695e7ef58932333fbeac582ae8bf107bc` |
| `src/baselines/multiris3d/private/b1_pai_build_initializations.m` | `0ecd2beed33a834ee6b147fc28a84f11578322209d63d988b896655ab474d34a` |
| `src/baselines/multiris3d/private/b1_pai_evaluate_candidates.m` | `33011bf5378125df4622850b22f84c1079d97c4e8c231bce197276fb041020c0` |
| `src/baselines/multiris3d/private/b1_pai_grid.m` | `ac834e9e55b0d4f534c1a8eb81931732f0a84e291ff5d987f4b292431001e3b7` |
| `src/baselines/multiris3d/private/b1_pai_objective.m` | `a6771bc801807c46e5949cfa8ada9c7cc082b8ace3c089dd39057a5fd73f1785` |
| `src/baselines/multiris3d/private/b1_pai_rank_initializations.m` | `898a71683cb93357c0fb2a86d98a5a9d9a6fe655a43b741da624df36bc6c0c16` |
| `src/baselines/multiris3d/estimate_channel_parameters.m` | `170e406a3ac281cc1cb146f139532ec9932dcba9edeb6339f37a1b9c3c1d68f0` |
| `src/baselines/multiris3d/estimate_delay_coarse.m` | `fb700248e6d072b1eebed28440bfdb917414431275acd62269cd53e8f7cb881b` |
| `src/baselines/multiris3d/estimate_spatial_coarse.m` | `6547dc0d27df3be5de32a05b7a6a2dd8651e6bbe7002c3eec54713557413317c` |
| `src/baselines/multiris3d/refine_channel_mle.m` | `d5b92a6e43130b48785f5616611e1b0a567cb2a90770242becbc953b32c5b540` |
| `src/baselines/multiris3d/estimate_position_coarse.m` | `176767693d77f8b802d225e8a94255c5495be0598979928235786817945e8cb6` |
| `src/baselines/multiris3d/refine_position_mle_eq57_identity.m` | `68b24daf57680a2eecc954963f57db4bc97f40f57534097149e0b585d8fa2a38` |
| `src/structured/current/signal/upa_steering.m` | `1fd2040e4862003a9fb156d8b61099b0305dd08501ea80abbc4e9982a7605fe4` |

已确认这些源码没有机器绝对路径、私有脚本调用、硬编码历史Trial/随机种子或真实位置答案。配置中的oracle搜索区域是既有实验信息条件，不将其表述为无真值辅助算法。
