# 修改清单

历史快照中的十个 `src/` 文件均未修改。复制前后的 SHA-256 相同。

新增文件：

- `experiment/paper_baseline_A_config.m`：核对历史固定物理配置，并执行已授权的 `TFT → TTT` 覆盖；定义真值附近搜索边界。
- `experiment/run_paper_baseline_A.m`：只使用 oracle support，调用原信道估计、粗定位和参数域单位权重精修。
- `experiment/separate_noiseless_observation.m`：调用原始无噪声物理生成器，并复用历史 DFT 分离步骤；不引入新物理模型。
- `experiment/static_check.m`：检查路径隔离、固定配置、1331 点网格和 Optimization Toolbox。
- `validation/run_all_validations.m`：实现获批的 A/B/C 验证、独立有限差分 Jacobian 检查、诊断字段和结果保存。
- `README.md` 与 `provenance/*`：运行说明和来源记录。

明确未做：

- 未修改或调用 VB bootstrap/后验更新。
- 未使用高斯位置先验、随机粗位置、运动轨迹或过程噪声。
- 未增加周期分支、UE 交换、多起点、失败回退或真值择优。
- 未约束连续精修搜索范围。
- 未新增比较2、理论界或正式 Monte Carlo。

环境提示：MATLAB R2022a 的 `restoredefaultpath` 报告安装目录内一个失效的 Add-Ons 路径及 locale 回退警告；项目函数解析、Optimization Toolbox 和全部验证不受影响。
