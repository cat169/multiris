# Baseline 2-S4独立正式接口与等价性审查

## 结论及版本边界

本次从已验证14候选分析实现中提取正式`infer_baseline2_frame_s4`，只执行S4所需计算。保留原MLE、Fisher、Fisher-MAP精修及候选选优规则；不调用Tracking临时分析脚本，不读取真实TX/RX位置，不在推断中写磁盘。

1600例留存审计通过；13个代表样本真实重算时，新S4与原14候选实现及历史检查点的四候选位置、完整目标、梯度、选择索引、迭代数、exitflag和异常标记一致。位置、目标和梯度最大差均为0。独立worktree仅使用本仓库三个生产路径，通过5次MAP交换样本检查；8项失败返回契约类测试通过。

审查分支`review/baseline2-s4-20261103`基于已审核B1-PAI提交`92cf12234d98ea3449568b15c86edaa654cf6c8c`。未修改B1-PAI、原B2-SR、原MLE/Fisher/MAP、Tracking正式入口或论文，未启动新的MC样本流、1600例全组重算或15,200例实验。

## 四候选：实际编号、初值、依赖和选择

令m=prior.mean，L=chol(prior.covariance,'lower')。内部编号严格保持原14候选MAT记录的前四项；历史CSV的candidateIndex为0起点，对应0/1/2/3。

| 内部编号/原名称 | 构造和优化 | 使用均值/TxX偏移 | 依赖 | 真值 |
|---|---|---|---|---|
| 1 originalSR | 原SR从m运行原MAP；若满足原交换触发条件，再从交换后的第一次终点运行原MAP，按原SR规则保留结果 | 首次从m；不加固定偏移 | 可能依赖第一次MAP终点做交换精修 | 不使用 |
| 2 mean | 独立重新从m执行原MAP | 是；无偏移 | 不依赖其他终点；与SR首次初值重复但按原路径保留重算 | 不使用 |
| 3 +TxX | 从m+L(:,1)独立执行原MAP | 是；4I下仅TxX+2m | 仅依赖固定先验/MLE/Fisher | 不使用 |
| 4 -TxX | 从m-L(:,1)独立执行原MAP | 是；4I下仅TxX−2m | 仅依赖固定先验/MLE/Fisher | 不使用 |

四个**输出候选**不是“四次优化永远独立”的同义词：候选1的SR交换保留必要依赖。候选2/3/4均可直接计算，不需要另外十个初值的结果。通常总计4次MAP；若SR交换触发，则5次。不能删除交换精修以强行保证4次，也不能将重复mean优化去重而改变原数值路径。1600例留存中1599例预计4次、1例5次；本轮Profiler实测12个代表例各4次、交换例5次。

“TxX”是原标签；任意相关协方差下L(:,1)可同时改变其他坐标，不能擅自改写成固定2m的单坐标扰动。1600例验证均为4I，本接口逐字保留实际L列构造，不据此宣称其他先验已有1600例验证。

正式主选择沿用原`selectedA(2)`：候选位置和报告目标必须有限、实数，独立复核目标与报告目标差≤1e-10*max(1,abs(reportedObjective))，按独立复核的完整J取最小；同值按1→4顺序取第一个。不用真实误差、不按exitflag或梯度删有限候选。原B规则（另要求正exitflag）的索引保留为`positiveExitflagIndex`诊断，不控制输出；1600例中A/B选择本来一致。合成边界测试确认有限exitflag=0候选仍可被主规则选中。

## 正式签名、输出和路径

```matlab
result = infer_baseline2_frame_s4( ...
    YSeparated, separatedNoiseVariance, baseReflections, ...
    positionPriorMessage, activeRisIndices, cfg);
```

与B2-SR相同六参数顺序。prior.mean为6×1列向量，prior.covariance为6×6正定Gaussian先验。调用方显式提供support和配置。本仓库生产路径：`baseline2/inference`、`src/baselines/multiris3d`、`src/structured/current/signal`。不需要Tracking目录；原MAP需要Optimization Toolbox的lsqnonlin。

主要输出：positionEstimate/staticPositionEstimate/mobilePositionEstimate、channel、fisher、map、objective、四项candidates、selectedIndex、positiveExitflagIndex、finite、converged、candidateStatus、exceptionIdentifier/Message、symmetryResolution、srOriginalException/srSwapException、callCounts和timing。

每候选保留label/initial/position、报告/独立完整目标及Fisher/先验分项、gradientInf、exitflag、iterations、finite、objectiveConsistent、异常、耗时和原map。`map.objective`保留原MAP报告值，`result.objective`为原A规则使用的独立复核值，正常允许浮点舍入差。finite、正exitflag、梯度收敛、真实误差分别评价，不混同。

## 完整依赖与计算流程

```
infer_baseline2_frame_s4
 → private/b2_s4_original_sr（从验证cached_sr提取；原SR数学和顺序不变）
   → estimate_channel_parameters（一次）
      → estimate_delay_coarse
      → estimate_spatial_coarse → upa_steering
      → refine_channel_mle → upa_steering
   → compute_channel_fim（每个active RIS一次，复用至所有MAP）
   → private/b2_s4_try_map → refine_position_map（首次及必要交换）
 → 原分析顺序的先验Cholesky/Fisher平方根
 → private/b2_s4_try_map → refine_position_map（mean/+TxX/-TxX三次）
 → private/b2_s4_residual（原multistart_residual仅更名提取，独立目标复核）
 → 主文件pick（原A规则最小完整目标）
```

私有残差是已验证独立复核函数的必要提取，没有重写实际优化目标；原`refine_position_map`内的目标、Jacobian和优化参数完全不改。原公开SR入口没有被本接口调用或修改；私有SR提取保留其实际计算顺序，新增的仅是耗时/计数与明确的异常返回处理。

本仓库已有upa_steering与历史Tracking版本仅换行字节不同，规范化后全文一致；代表重算和独立checkout均使用本仓库版本，MLE/Fisher逐位一致。

## 冻结的Fisher-MAP数学

MLE使用原时延粗搜、0.02空间网格及single-start fminsearch连续精修，参数不变。每个可见RIS由估计q=[xi,zeta,Delta_f*tau]、rho、W和分离后复噪声方差计算

`F=2*abs(rho)^2/nv*real(D'*D-(D'*mu)*(D'*mu)'/(mu'*mu))`。

未知复增益按原Schur消元，不显式构造大投影阵，不求F逆/伪逆、不加ridge。原对称化/eig规则保留：超过100*eps*max(norm(F,2),realmin)的负特征值报错，舍入范围内负值截0，R=diag(sqrt(e))*V'。

`J(x)=0.5*sum_l r_l(x)'*F_l*r_l(x)+0.5*(x-m)'*C^(-1)*(x-m)`。

级联几何及周期[2;2;1]不变；每RIS在最近分支周围枚举{-1,0,1}^3，取Fisher加权代价最小的残差。先验Cholesky三角解、解析Jacobian不变。原lsqnonlin使用trust-region-reflective，MaxIterations=500、MaxFunctionEvaluations=5000、FunctionTolerance=1e-12、StepTolerance=1e-12、OptimalityTolerance=1e-10，无位置边界、额外惩罚或真值约束。

## 正式失败返回契约

- 四个均满足原A规则：candidateStatus=allCandidatesValid。
- 部分无效但仍有候选：partialCandidateFailure，按同一规则选择有效候选，完整保留失败候选证据。
- 全部无效：allCandidatesFailed，finite=false，selectedIndex=0，位置/目标NaN、map为空struct，异常baseline2S4:NoValidCandidate；不伪造正常结果。
- 原MAP已明确的几何域异常`refine_position_map:PositionAtRisCenter`可记为候选失败后继续。原SR首次遇到它时仍运行另外三初值；依赖交换精修遇到它时保留原SR首次结果，并单独记录srSwapException。
- 错误先验、非PSD Fisher、未知MATLAB/工具箱/程序异常继续抛出，不统一吞掉。此保守边界不声称识别所有可能的工具箱数值错误。

这是本轮明确新增的失败契约：历史14候选脚本对额外初值捕获所有异常，而SR异常会直接终止验证。新接口区分已知几何失败与程序错误；1600例原无候选异常，成功数值路径和原选优规则不受影响。8项类测试使用隔离solver doubles验证上述返回，不将其冒充数值回归。

## 1600例留存结果审计（不是1600例重算）

实际加载旧checkpoint.mat的records/completed/spec，1600/1600完成。逐Trial核对前四候选label、六维初值、有限性和异常标记；从四个留存完整目标重新取最小，选中索引与原selectedA(2)一致，并重建最终六维位置；B规则索引也与原selectedB(2)一致。原检查点、CSV和1600例验证源码均保持原SHA-256。

| 场景 | SNR | 留存例数 |
|---|---:|---:|
| 011 | 0 | 300 |
| 011 | 10 | 300 |
| 011 | 20 | 300 |
| 111 | 10 | 100 |
| 111 | 20 | 100 |
| 110 | 10 | 100 |
| 110 | 20 | 100 |
| 101 | 10 | 100 |
| 101 | 20 | 100 |
| 111 | -15 | 100 |

留存审计输出`test_results/baseline2_s4_interface/retained_1600_audit.csv`保留1600例最终位置、四候选目标、exitflag及选中编号；原始完整候选终点仍在未修改的旧检查点。

## 13例真实重算与固定容差

重算覆盖全部10个已验证场景/SNR单元，并额外覆盖011/10dB Trial27/209和111/-15dB Trial81。原1600检查点未保存Y，因此按原配置、幅度、rawNoiseVariance和原signalSeed确定性重建相同观测；重新MLE与留存channel逐位一致，Fisher也逐位一致。原14候选函数和新接口收到同一份Y/nv/W/prior。没有生成新的随机样本或重跑VB/Baseline1/PEB。

事先固定：初值及Fisher精确相等；位置无穷范数差≤1e-8m；完整目标差≤1e-10*max(1,abs(reference))；梯度差≤1e-8*max(1,abs(reference))；候选索引、exitflag、迭代数、异常标记精确相等。实际四候选位置、目标和梯度差最大值均0，最终位置与原分析逐位一致。

| 场景 | SNR | Trial | 选中内部编号 | MAP次数 | 最大位置/目标差 | 最终Joint误差(m，仅评价) |
|---|---:|---:|---:|---:|---:|---:|
| 011 | 0 | 1 | 3 | 4 | 0 / 0 | 0.433119839047 |
| 011 | 10 | 1 | 4 | 4 | 0 / 0 | 0.361793747877 |
| 011 | 20 | 1 | 3 | 4 | 0 / 0 | 0.061892309051 |
| 111 | 10 | 1 | 4 | 4 | 0 / 0 | 0.015592495545 |
| 111 | 20 | 1 | 3 | 4 | 0 / 0 | 0.008300861469 |
| 110 | 10 | 1 | 3 | 4 | 0 / 0 | 0.031084153670 |
| 110 | 20 | 1 | 4 | 4 | 0 / 0 | 0.022087749456 |
| 101 | 10 | 1 | 4 | 4 | 0 / 0 | 0.025461876491 |
| 101 | 20 | 1 | 4 | 4 | 0 / 0 | 0.008050535791 |
| 111 | -15 | 1 | 4 | 4 | 0 / 0 | 0.583164029195 |
| 011 | 10 | 27 | 3 | 4 | 0 / 0 | 0.061556842564 |
| 011 | 10 | 209 | 3 | 4 | 0 / 0 | 0.082682230584 |
| 111 | -15 | 81 | 4 | 5 | 0 / 0 | 45.474589898100 |

Profiler实际观察到13次MLE、各样本每active RIS一次Fisher，新S4共53次MAP（12×4+1×5），没有调用其余十个固定初值。接口返回计数与Profiler核对。独立worktree采用原Trial81冻结输入，production路径只含本仓库，检查到11个用户源码依赖、无Tracking依赖并得到相同终点/选中编号及5次MAP。8项类测试在独立worktree也全部通过。

测试脚本：`tests/baseline2_s4/run_baseline2_s4_equivalence.m`、`check_baseline2_s4_standalone.m`、`Baseline2S4ContractTest.m`和隔离fixture。必要大型检查点/观测不提交；读取历史外部验证函数仅发生在等价性脚本，正式接口没有此依赖。MATLAB R2022a；主文件及等价性脚本静态检查无问题。本机restoredefaultpath的既有addons/Perl locale告警没有阻断核验。

## 源码SHA-256与提交集合

仅新增四个生产文件、更新本文、加入轻量验证脚本及fixture。七个原计算依赖已在基准仓库中，不重复提交。新文件采用LF；已有文件SHA按当前工作区原始字节列出，Git checkout换行规范化可改变字节SHA，全文和数学内容仍应一致。

| 文件 | SHA-256（本地原始字节） | 提交状态 |
|---|---|---|
| `baseline2/inference/infer_baseline2_frame_s4.m` | `95dc091a91cd67e8fb5b8d8df64fe1b4aadab3646392b9290d64f066ece00fc6` | 新增 |
| `baseline2/inference/private/b2_s4_original_sr.m` | `20c5420027ebbc20285cfe241d0d12ea8ff0c1419ff72f54ad328f1b6c28c96f` | 新增 |
| `baseline2/inference/private/b2_s4_try_map.m` | `e99459ae9ccb7c13b872a78debe2b0d8a91051c4ddeb1db70d21f1a075f06006` | 新增 |
| `baseline2/inference/private/b2_s4_residual.m` | `91175b6b61ca126b81912fe2531457d36a4e8cf9dc648d2dd5d4d31347cdb156` | 新增 |
| `baseline2/inference/compute_channel_fim.m` | `9262670d5461c1a1ef65814ee5cddfe8a5de38ff2f03cda69f87ad7accef761d` | 基准已有、未修改 |
| `baseline2/inference/refine_position_map.m` | `b5fb54fd390195cd9587762368ba318518a12d9cc9e2662465b66c0830c0c246` | 基准已有、未修改 |
| `src/baselines/multiris3d/estimate_channel_parameters.m` | `170e406a3ac281cc1cb146f139532ec9932dcba9edeb6339f37a1b9c3c1d68f0` | 基准已有、未修改 |
| `src/baselines/multiris3d/estimate_delay_coarse.m` | `fb700248e6d072b1eebed28440bfdb917414431275acd62269cd53e8f7cb881b` | 基准已有、未修改 |
| `src/baselines/multiris3d/estimate_spatial_coarse.m` | `6547dc0d27df3be5de32a05b7a6a2dd8651e6bbe7002c3eec54713557413317c` | 基准已有、未修改 |
| `src/baselines/multiris3d/refine_channel_mle.m` | `d5b92a6e43130b48785f5616611e1b0a567cb2a90770242becbc953b32c5b540` | 基准已有、未修改 |
| `src/structured/current/signal/upa_steering.m` | `1fd2040e4862003a9fb156d8b61099b0305dd08501ea80abbc4e9982a7605fe4` | 基准已有、未修改 |

本轮另外核验218个受保护文件SHA未变，含两工程原源码、旧分析函数和历史独立验证数据；两个原工作树HEAD和暂存区也未变化。原未提交MAT修改、历史删除及B1测试均保留。完整新MAT诊断只保存在本地test_results，不进入Git。

## 剩余问题与适用边界

1. S4不是全局最优保证。Trial81低SNR错解仍为45.474589898100m，更低MAP目标并不保证真实误差更小，没有以真值强行拒绝它。
2. 位置有限、正exitflag与梯度达标不同；目标/步长停止不保证严格梯度收敛。原SR内部的正exitflag交换门限保留，与S4外层A规则区分。
3. 重复mean求解、每次MAP重做Fisher eig/先验Cholesky和独立目标复核均保留，未做会改变执行路径的优化。推断本身没有保存文件或大量网格诊断。
4. 四候选中的SR可能需要依赖性第五次MAP；这不是提取错误，而是原策略行为。
5. 新失败返回边界已经类测试，但原1600例没有候选异常，不能将留存样本称为覆盖全部真实数值故障。
6. 原S4方向选择来自历史异常，存在事后方向选择偏差；本次工程等价性不增加其统计泛化证据，也不等于全部1600例重新运行新接口。

截至本审查包，没有改动Tracking正式入口；等待源码审查后再单独开展集成。
