function [YSeparated, separatedNoiseVariance, baseReflections] = ...
    generate_noisy_separated_signal(cfg)
%GENERATE_NOISY_SEPARATED_SIGNAL 生成并分离一帧多 RIS 加噪信号。
%{
输入 cfg 为场景、RIS、OFDM、路径可见性和 SNR 配置。函数先调用
generate_noiseless_signal 生成原始无噪信号与 DFT 码本，再按显式信号功率
加入圆对称复高斯白噪声，最后用 DFT 码字投影得到逐 RIS 加噪观测。

输出 YSeparated 的维度为“子载波 × 基础配置 × RIS”；
separatedNoiseVariance 是投影后每个复样本的总噪声方差；baseReflections
保留生成本帧观测所使用的各 RIS 基础反射配置。
%}

[Y, ~, risReflections, codebook] = generate_noiseless_signal(cfg);

signalPower = mean(abs(Y(:)).^2);
if signalPower == 0
    noisySignal = Y;
    noiseVariance = 0;
else
    noiseVariance = signalPower / ...
        10^(cfg.noise.snr_db / 10);
    noise = sqrt(noiseVariance / 2) * ( ...
        randn(size(Y)) + 1i * randn(size(Y)));
    noisySignal = Y + noise;
end

% 按 DFT 块重排加噪信号，并投影到每条 RIS 的正交码字。
numSubcarriers = cfg.ofdm.num_subcarriers;
baseConfigCount = cfg.dft.base_configs_per_block;
numBlocks = cfg.dft.num_blocks;
numRis = cfg.ris.count;

baseReflections = complex(zeros( ...
    cfg.ris.num_elements, baseConfigCount, numRis));
for ell = 1:numRis
    baseReflections(:, :, ell) = ...
        risReflections(:, 1:baseConfigCount, ell) ./ ...
        conj(codebook(1, ell));
end

signalBlocks = reshape( ...
    noisySignal, numSubcarriers, baseConfigCount, numBlocks);
YSeparated = complex(zeros( ...
    numSubcarriers, baseConfigCount, numRis));

for ell = 1:numRis
    for block = 1:numBlocks
        YSeparated(:, :, ell) = ...
            YSeparated(:, :, ell) + ...
            codebook(block, ell) * ...
            signalBlocks(:, :, block) / numBlocks;
    end
end

% 独立 DFT 块经归一化相加后，总复噪声方差缩小 numBlocks 倍。
separatedNoiseVariance = noiseVariance / numBlocks;
end
