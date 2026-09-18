# 静态检查报告

检查日期：2026-09-19。

本检查未启动 MATLAB、Octave、验证、仿真或实验监控。

已检查：

- 三个 `.m` 文件均以函数声明开始、以 `end` 结束。
- 圆括号和方括号文本计数平衡。
- Git whitespace 检查无错误。
- 固定 SNR 为 `-20:5:30`，固定 seed 为 `20261001:20261300`。
- 正式入口在每个 SNR/trial 前使用固定 seed 重置 RNG。
- 检查点按 SNR 独立命名，并校验 schema、完整配置、源码哈希、SNR、trial 数、seed 和数组尺寸。
- 汇总与绘图文件不包含信号生成、信道估计、定位或 VB 调用。
- MATLAB 文件中未发现 bootstrap、后验、先验惩罚、交换、多起点、失败回退、周期分支、理论界或 BCRB 增强。
- 正式结果写入 `formal_results`，不会覆盖已有 `results` 验证目录。
- 正式入口使用 `onCleanup` 恢复调用前的 MATLAB path 和 RNG，不调用 `savepath`。
- 目标目录在写入前不存在 `formal` 和 `formal_results`，没有覆盖同名用户文件。

限制：按照用户要求，本轮没有使用 MATLAB `checkcode` 或实际函数解析。因此 MATLAB 语法、Java SHA-256 调用、检查点写入及绘图导出仍属于“尚未执行验证”。
