# B1-PAI 最小工程修复与回归验证

## 结论

本次审查指定的候选异常处理、外部MLE结构检查、Tx零距离和重复rank问题已修复。历史162例失败和100例正常回归全部通过；最终位置、所有候选位置和初/终目标保持一致，选中候选、exitflag、迭代数、函数调用数、MLE和质量标签不变。历史失败样本Joint RMSE为 **61.3577904760193 m**。

**可进入下一轮Tracking正式集成，但本轮没有修改Tracking入口。** 集成调用方必须适配空`coarse/refined`、`finite=false`及返回异常标识，并继续独立捕获程序/配置错误。此工程结论不代表先前已报告的低SNR远端解、切向方向或RIS中心数值病态得到算法性改善。

## 1. 修改范围与原因

基于 `fe8253e8f9885af896cf7c037b2fc86427194c14`，生产代码只修改：

| 文件 | 发现的问题 | 最小修复 |
| --- | --- | --- |
| `infer_baseline1_frame_projection_assisted.m` | 外部MLE可带错误支持、尺寸、非有限值或不匹配噪声 | 增加局部检查函数，仅非空第5参数调用时执行；四参数MLE路径不变 |
| `private/b1_pai_evaluate_candidates.m` | NaN/Inf目标可触发一致性断言并终止整帧；旧catch可吞程序错误 | 先检查有限实数性；无效数值原样记录并继续；一致性断言容差不变；catch只包围原精修调用且仅白名单数值异常可恢复 |
| `private/b1_pai_build_initializations.m` | Tx恰在RIS中心会除零；同一pair重复rank | 仅距离精确等于0时跳过Tx网格、记录`txAtRisCenter`；复用第一次rank值 |

原信道MLE、原粗定位、原精修、网格生成、原目标表达式和排序/去重helper均未修改。新增测试位于 `tests/baseline1_pai`，不插入生产推断路径。更新源码审查文档，保留原有实验文件及未提交修改。

## 2. 关键diff与不变量

候选处理改为：

```matlab
try
    f = refine_position_mle_eq57_identity(...); % 原调用
catch ex
    if ~any(strcmp(ex.identifier,numericalIds)), rethrow(ex); end
    % 只记录明确的单候选数值错误，继续下一候选
end
% 原目标表达式求值后：
if ~q.valid
    q.exceptionIdentifier = 'baseline1PAI:NonfiniteCandidate';
    % 保留NaN/Inf/非实数原值及具体不合格字段，不执行一致性断言
    continue;
end
% 原Initial/FinalObjectiveMismatch断言及1e-12*max(1,J)容差不变
```

恢复白名单仅包含本机原求解路径可识别的：

- `refine_position_mle_eq57_identity:PositionAtRisCenter`
- `optimlib:commonMsgs:UsrObjUndefAtX0`
- `optimlib:commonMsgs:UserJacUndefAtX0`

后两项来自本机R2022a `lsqncommon` 的初始残差/Jacobian非有限检查。其他异常（包括维度、类型、配置和有限目标不一致）不吞掉。数值标识及消息写入候选记录；没有按消息文本模糊分类。全部无效时保持原顶层失败结构，不向上抛`NoValidRefinement`。

外部MLE必须为标量结构并有4个必要字段；支持集合及顺序与配置一致，估计是长度匹配的有限实数向量，存在噪声字段时与本次输入完全相同。不重新执行MLE，**这些检查不能证明与观测严格同源**，仍须调用方保证。没有增设输入真值、先验、Fisher或误差选择规则。

几何改动为 `distance=norm(d); if distance==0, ... continue; end`，不使用新非零距离门限。pair仍按`rank<2`拒绝，只复用`diag.pairRank`。投影公式、局部轴、正法向、RIS配对、line-i最近点、等权融合、原目标排序、1e-6m去重及最多8初值均冻结。有限exitflag<=0或大梯度候选不自动剔除，质量标签不参与最小原目标选择。

## 3. 逐Trial与逐候选回归

复用既有验证checkpoint和历史方案A缓存。前202例直接复用原观测、反射矩阵和原MLE；后60例按原种子重建观测，仍复用该Trial原MLE。没有新种子Monte Carlo、VB、B2、PEB或MLE重推断。通过公共推断接口调用当前生产文件，不使用替身。

| 路径 | 历史失败 | 历史正常 | 中高SNR正常 | 最大最终位置差/m | 最大候选位置/目标差 | 编号/停止/质量变化 |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| 111 | 96/96 | 10/10 | 15/15 | 0 | 0 | 0 |
| 110 | 33/33 | 10/10 | 15/15 | 0 | 0 | 0 |
| 101 | 10/10 | 10/10 | 15/15 | 0 | 0 | 0 |
| 011 | 23/23 | 10/10 | 15/15 | 0 | 0 | 0 |

中高SNR为各状态0/10/20dB各5个原正常Trial。正常100例均不投影、额外精修0次、原几何输出一致；原正常分支仍仅1次精修。默认四参数一次MLE、一次粗定位、一次精修另由隔离计数测试确认。没有用RMSE一致代替逐Trial检查。

每个候选均比较六维初/终点、初/终完整目标、候选可用性、选中编号、exitflag、迭代次数、函数调用次数、一阶指标、最近RIS距离及完整质量原因结构。沿用位置绝对1e-8m、目标1e-12*max(1,|J|)、一阶指标1e-8*max(1,|g|)、RIS距离1e-8m容差；本轮实际上述数值差均为0。逐Trial262行结果及当前候选保存在本地 `test_results/baseline1_pai_hardening/regression.csv`、`regression_candidates.mat`，不随源码提交。

## 4. 隔离边界测试

23项测试在MATLAB MCP中 **23 Passed / 0 Failed / 0 Incomplete**；完整回归runner末尾再次执行，同样23/23通过。`matlab.unittest.TestCase`只测试公共接口。fixture把生产源码复制到临时隔离目录，在该副本使用精修/粗定位/MLE替身；生产原精修没有被修改。每项结束恢复路径并清理自己创建的目录。

| Test | Result |
| --- | --- |
| `rejectNonstructChannel` | Passed |
| `rejectMissingField` | Passed |
| `rejectSupportOrder` | Passed |
| `rejectWrongVectorLength` | Passed |
| `rejectMatrixEstimate` | Passed |
| `rejectNaN` | Passed |
| `rejectInf` | Passed |
| `rejectComplex` | Passed |
| `rejectNoiseMismatch` | Passed |
| `skipExactRisCenter` | Passed |
| `continueAfterNaNPosition` | Passed |
| `continueAfterInfFinalObjective` | Passed |
| `continueAfterInfInitialObjective` | Passed |
| `allCandidatesInvalidReturnsFailure` | Passed |
| `allGridPointsInvalidReturnsFailure` | Passed |
| `initialMismatchIsFatal` | Passed |
| `finalMismatchIsFatal` | Passed |
| `dimensionErrorIsFatal` | Passed |
| `unknownErrorIsFatal` | Passed |
| `knownNumericalExceptionContinues` | Passed |
| `nonconvergenceDoesNotFilterFiniteCandidate` | Passed |
| `fourArgumentPathEstimatesMleOnce` | Passed |
| `threeReturnCasesHaveSameTopLevelFields` | Passed |

三类语义明确分开：错误外部输入预期拒绝；NaN/Inf或白名单数值异常是可继续的单候选失败；有限目标不一致、维度及未知程序异常必须立即报告。全部候选无效的返回断言覆盖`fallbackTriggered=true`、`finite=false`、NaN位置、`selectedIndex=0`、`baseline1PAI:NoValidRefinement`及空`coarse/refined`。正常/回退成功/全部失败三者顶层字段集合一致。

## 5. MATLAB路径、语法和依赖

本机MATLAB R2022a。生产入口、原信道估计、原粗定位及原精修解析到本仓库 `src/baselines/multiris3d`，未解析到测试替身。三个修改文件的 `checkcode(...,'-id')`结果均为空（0问题）；`matlab.codetools.requiredFilesAndProducts`解析成功，依赖记录保存在本地测试结果目录。测试替身只在边界测试临时镜像可见。

原MLE、几何、精修和signal steering没有修改；181个既有受保护源码/历史数据/先前验证文件SHA-256与已有清单一致。Tracking代码与历史报告均只读，未写入或覆盖。

## 6. 修复前后源码SHA-256

修复前为初版审查提交的原始Git blob字节；修复后为当前验证工作文件字节，按原字节暂存。换行不作为算法差异，未修改的3个helper字节哈希保持相同。

| File | Before SHA-256 | After SHA-256 |
| --- | --- | --- |
| `src/baselines/multiris3d/infer_baseline1_frame_projection_assisted.m` | `084040ad2c79683ef6eddf520547212e9770489d881629e4c8b966ac92831411` | `06fe32fc43d79c3f4ca4192e529b1a8695e7ef58932333fbeac582ae8bf107bc` |
| `src/baselines/multiris3d/private/b1_pai_evaluate_candidates.m` | `dd0d43b38f99e9db21bf9919b73ec5032c1c1701a98f4cf912970fa9e0364178` | `33011bf5378125df4622850b22f84c1079d97c4e8c231bce197276fb041020c0` |
| `src/baselines/multiris3d/private/b1_pai_build_initializations.m` | `42a9fe1467726217693ce0398c1b7f9747b4efd0f4bfe13a37a7a53c7329bfc9` | `0ecd2beed33a834ee6b147fc28a84f11578322209d63d988b896655ab474d34a` |
| `src/baselines/multiris3d/private/b1_pai_grid.m` | `ac834e9e55b0d4f534c1a8eb81931732f0a84e291ff5d987f4b292431001e3b7` | `ac834e9e55b0d4f534c1a8eb81931732f0a84e291ff5d987f4b292431001e3b7` |
| `src/baselines/multiris3d/private/b1_pai_objective.m` | `a6771bc801807c46e5949cfa8ada9c7cc082b8ace3c089dd39057a5fd73f1785` | `a6771bc801807c46e5949cfa8ada9c7cc082b8ace3c089dd39057a5fd73f1785` |
| `src/baselines/multiris3d/private/b1_pai_rank_initializations.m` | `898a71683cb93357c0fb2a86d98a5a9d9a6fe655a43b741da624df36bc6c0c16` | `898a71683cb93357c0fb2a86d98a5a9d9a6fe655a43b741da624df36bc6c0c16` |

## 7. 是否需要重跑800例

**本轮不需要重跑800例独立优化，有核查依据而非默认免检。** 已只读审计800条留存记录及961个候选：外部MLE支持/长度/有限性及噪声全部通过；原G0的Tx网格均不与有效RIS中心重合；没有非有限候选，也没有旧版吞掉而新版会转为fatal的异常。新增守卫对该留存集合不改变计算分支，rank复用不改判据；再结合262例逐候选真实重算和边界替身测试，已覆盖本次工程变更。

只读审计不等于重新运行800例，也不证明未来任意配置永无数值异常。如后续改变求解器、几何、初值或白名单，应按变更重新评估。其结果保存在本地 `frozen_800_guard_audit.csv`；本轮未重新生成800例观测或进行其位置优化。

## 8. 剩余契约与集成结论

- 本次指定工程问题已消除，有限输出逐Trial及逐候选与原验证版一致，可进入下一轮正式集成。
- 观测同源性仍由调用方保证；沿用原配置schema，包括现有`geometry.static_user/mobile_user`字段。新检查不声称验证了外部MLE的统计正确性。
- 未知/其他MATLAB版本的数值异常保守地继续抛出，不保证所有异常都可恢复，也不应按错误文本自动吞掉程序错误。外层Monte Carlo仍需独立异常处理。
- 空结构、有限性、原异常、回退状态和质量标签必须在Tracking适配中保存；不能无条件读取空`refined.exitflag`。
- 原算法的低SNR误差、切向方向、负射线与RIS中心邻域结果没有被“修复”；仍须报告有限率、质量与全Trial分母定位成功率。

## 9. 可重复运行命令

将本仓库根目录设为MATLAB当前目录，Tracking为同级目录时：

```matlab
addpath(fullfile(pwd,'tests','baseline1_pai'));
results = runtests(fullfile('tests','baseline1_pai','Baseline1PaiHardeningTest.m'));
assert(all([results.Passed]));
trackingRoot = fullfile(fileparts(pwd),'Tracking_code');
outputDirectory = fullfile(pwd,'test_results','baseline1_pai_hardening');
run_b1_pai_hardening_validation(trackingRoot,outputDirectory);
audit_frozen_independent_guards(trackingRoot,outputDirectory);
```

只有回归需要既有本地MAT输入；边界测试自包含。新测试结果不覆盖历史实验。源码/必要测试/文档推送至原审查分支，不合并main，不修改Tracking正式入口。
