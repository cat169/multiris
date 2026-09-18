function cfg = simulation_config()
%SIMULATION_CONFIG 配置单帧推断和多帧真实数据生成。
%{
集中定义场景几何、RIS/OFDM 参数、仿真真值、初始化先验、局部 VB 停止阈值、
噪声和多帧随机游走设置。
%}

%% 1. 物理常数

cfg.constants.c = 3e8;
cfg.constants.lambda = 10.7e-3;
cfg.constants.fc = ...
    cfg.constants.c / cfg.constants.lambda;


%% 2. 三维场景参考范围

cfg.geometry.room_min = [-15; -15; 0];
cfg.geometry.room_max = [ 15;  15; 6];

%% 3. RIS 阵列参数

cfg.ris.count = 3;

% nx 对应物理局部水平阵列维度 Ny；
% ny 对应物理局部竖直阵列维度 Nz。
cfg.ris.nx = 10;
cfg.ris.ny = 10;
cfg.ris.num_elements = ...
    cfg.ris.nx * cfg.ris.ny;

% 半波长阵元间距。
cfg.ris.dx = cfg.constants.lambda / 2;
cfg.ris.dy = cfg.constants.lambda / 2;


%% 4. 各 RIS 的位置

% 每一列对应一块 RIS。
cfg.geometry.ris_centers = [
     0,   0, -15;
    20, -20, 5;
    10,  10,  4.5
];



%% 5. RIS 的局部阵列轴

% 第一局部阵列轴控制组合空间频率 xi：
cfg.geometry.ris_axis_x = [
     1, -1, 0;
     0,  0, -1;
     0,  0, 0
];

% 第二局部阵列轴控制组合空间频率 zeta：
cfg.geometry.ris_axis_y = [
     0, 0, 0;
     0, 0, 0;
     1, 1, 1
];

%% 6. RIS 的正法向轴

% 每一列对应一块 RIS。
% 正法向定义为 RIS 面向有效定位区域的一侧。
%
% RIS 1: [0; -1; 0]
% RIS 2: [0;  1; 0]
% RIS 3: [1;  0; 0]
%
% 当前 static/mobile UE 均位于三块 RIS 的正法向半空间，
% 因此后续由 xi、zeta 恢复第三个局部方向分量时统一取正根。

cfg.geometry.ris_axis_normal = [
     0,  0,  1;
    -1,  1,  0;
     0,  0,  0
];
%% 6. 仿真真值位置

cfg.geometry.static_user = [-8; -6; 1.5];
cfg.geometry.mobile_user = [ 7;  7; 1.5];
%% 7. OFDM 参数

cfg.ofdm.num_subcarriers = 40;
cfg.ofdm.bandwidth = 100e6;

cfg.ofdm.subcarrier_spacing = ...
    cfg.ofdm.bandwidth / ...
    cfg.ofdm.num_subcarriers;


%% 8. RIS 的 DFT 时间编码

% 三块 RIS 使用长度为 3 的正交 DFT 编码。
cfg.dft.num_blocks = cfg.ris.count;

% 每个 DFT 块包含 40 个基础 RIS 配置。
cfg.dft.base_configs_per_block = 40;

% 总 OFDM 符号数为 3 × 40 = 120。
cfg.dft.num_symbols = ...
    cfg.dft.num_blocks * ...
    cfg.dft.base_configs_per_block;


%% 9. 第一帧位置初始化

% 两个用户的位置都未知，采用相同的各向同性初始化误差。
cfg.initialization.position_variance = 0.8;

cfg.initialization.static_position_covariance = ...
    cfg.initialization.position_variance * eye(3);

cfg.initialization.mobile_position_covariance = ...
    cfg.initialization.position_variance * eye(3);


%% 10. 路径可见性

% 真值。
% cfg.paths.visible = true(cfg.ris.count, 1);
cfg.paths.visible = [true;false;true];
% 多帧真实路径与推断侧帧间预测使用相同的二状态 Markov 参数。
cfg.paths.birth_probability = 0.9;
cfg.paths.death_probability = 0.05;
cfg.paths.decision_threshold = 0.5;
cfg.paths.markov_enabled = true;
% 第一帧采用无信息 Bernoulli 预测先验 pi^-=0.5。
cfg.local_vb.predicted_visibility_probability = ...
    0.5 * ones(cfg.ris.count, 1);


%% 11. 路径增益先验

% 推断固定使用方差为 1e-17 的宽零均值复高斯先验：
%
% p(rho_l) = CN(0, v_rho0_l)。
%
% 仿真真值 rho 不从该先验抽样，而是在 generate_noiseless_signal
% 中根据静态用户—RIS—移动用户的距离公式确定。
cfg.local_vb.path_gain_prior_variance = [ ...
    1e-17; ...
    1e-17;1e-17];


%% 12. 单帧局部 VB 迭代设置

% 每块 RIS 独立执行局部坐标上升 VB。
cfg.local_vb.max_iterations = 50;

% 使用 VM 第一圆周矩的变化判断三个相位是否收敛。
cfg.local_vb.circular_moment_tolerance = 1e-6;

% 路径增益后验均值和方差的相对变化阈值。
cfg.local_vb.gain_relative_tolerance = 1e-6;

% 路径存在软后验概率 E[z] 的轮间绝对变化阈值。
cfg.local_vb.visibility_probability_tolerance = 1e-6;

%% 13. 噪声设置

cfg.noise.snr_db = 10;

%% 14. 多帧真实轨迹生成

% 数组第 1 帧对应数学时刻 t=0；单位 m^2。
cfg.tracking.num_frames = 300;
cfg.tracking.mobile_process_covariance = 0.05 * eye(3);

end
