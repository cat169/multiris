# B1-PAI 工程修复验证

在仓库根目录作为 MATLAB 当前目录时可执行：

```matlab
addpath(fullfile(pwd,'tests','baseline1_pai'));
results = runtests(fullfile('tests','baseline1_pai','Baseline1PaiHardeningTest.m'));
assert(all([results.Passed]));
trackingRoot = fullfile(fileparts(pwd),'Tracking_code');
outputDirectory = fullfile(pwd,'test_results','baseline1_pai_hardening');
run_b1_pai_hardening_validation(trackingRoot,outputDirectory);
audit_frozen_independent_guards(trackingRoot,outputDirectory);
```

边界测试自身不需要历史数据。`PaiHardeningFixture` 将当前源码复制到独立临时目录，在复制环境用 `.fixture` 文件替代精修、粗定位和MLE；只通过公共推断接口测试，不改生产精修。测试结束恢复路径并删除自己创建的临时目录。

回归需要已保存的 Tracking `reports/baseline1_pai_validation/checkpoint.mat` 与 `reports/baseline1_fallback_comparison/checkpoint.mat`，不随源码提交。262例全部复用原MLE；前202例复用原观测，后60例按原种子重建观测、仍复用原MLE。没有新种子仿真、VB、B2或800例重新优化。新输出位于本仓库 `test_results`，不覆盖历史结果，也不修改 Tracking 文件。

回归逐候选比较初/终目标、位置、选中编号、exitflag、迭代数、函数调用数、质量标记和MLE。预先沿用绝对位置1e-8m、目标1e-12*max(1,|J|)、一阶指标1e-8*max(1,|g|)、RIS距离1e-8m；整数状态/质量/初始化完全一致。RMSE仅额外核验，不代替逐Trial核对。

冻结800例检查只审计新增分支谓词，不重新求解：外部MLE、噪声匹配、Tx网格零距离、非有限候选及旧版吞掉但新版会传播的异常。若任一检查失败，应检查受影响样本，不能据此直接宣称无需回归。
