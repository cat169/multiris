function [Y, path, risReflections, codebook] = ...
    generate_noiseless_signal(cfg)
%GENERATE_NOISELESS_SIGNAL 生成一帧多 RIS 原始无噪信号。
%{
输入 cfg 为场景、RIS、OFDM 和路径可见性配置。函数依次计算路径参数、生成
DFT 编码反射配置，并叠加各 RIS 路径得到原始无噪信号。

输出 Y 为“子载波 × OFDM 符号”的原始无噪信号；path、risReflections 和
codebook 分别保留路径参数、完整反射配置和 DFT 码本。
%}

% 根据两段用户—RIS 几何计算级联距离、时延、局部空间参数和复路径增益。
deltaStatic = cfg.geometry.static_user - cfg.geometry.ris_centers;
deltaMobile = cfg.geometry.mobile_user - cfg.geometry.ris_centers;

dStatic = vecnorm(deltaStatic, 2, 1);
dMobile = vecnorm(deltaMobile, 2, 1);
unitStatic = deltaStatic ./ dStatic;
unitMobile = deltaMobile ./ dMobile;

distance = dStatic + dMobile;
directionSum = unitStatic + unitMobile;

path.d_static = dStatic.';
path.d_mobile = dMobile.';
path.distance = distance.';
path.tau = path.distance / cfg.constants.c;
path.xi = sum( ...
    cfg.geometry.ris_axis_x .* directionSum, 1).';
path.zeta = sum( ...
    cfg.geometry.ris_axis_y .* directionSum, 1).';
path.rho = ( ...
    cfg.constants.lambda^2 ./ ...
    (16 * pi^2 * dStatic .* dMobile) .* ...
    exp(-1i * 2 * pi / cfg.constants.lambda * distance)).';

% 生成基础随机单位模反射配置，并用正交 DFT 码字扩展到全部时隙。
numElements = cfg.ris.num_elements;
numRis = cfg.ris.count;
numBlocks = cfg.dft.num_blocks;
baseConfigCount = cfg.dft.base_configs_per_block;
numSymbols = cfg.dft.num_symbols;

blockIndex = (0:(numBlocks - 1)).';
risIndex = 0:(numRis - 1);
codebook = exp( ...
    -1i * 2 * pi / numBlocks * (blockIndex * risIndex));

baseReflections = exp( ...
    1i * 2 * pi * rand(numElements, baseConfigCount, numRis));

risReflections = complex(zeros(numElements, numSymbols, numRis));
for ell = 1:numRis
    for block = 1:numBlocks
        columns = (block - 1) * baseConfigCount + (1:baseConfigCount);
        risReflections(:, columns, ell) = ...
            conj(codebook(block, ell)) * ...
            baseReflections(:, :, ell);
    end
end

% 叠加各 RIS 的时延响应、UPA 空间响应、路径可见性和复增益。
numSubcarriers = cfg.ofdm.num_subcarriers;
subcarrierIndex = (0:(numSubcarriers - 1)).';
elementIndexX = (0:(cfg.ris.nx - 1)).';
elementIndexY = (0:(cfg.ris.ny - 1)).';
Y = complex(zeros(numSubcarriers, numSymbols));

for ell = 1:numRis
    delayResponse = exp( ...
        -1i * 2 * pi * subcarrierIndex * ...
        cfg.ofdm.subcarrier_spacing * path.tau(ell));

    steeringX = exp(1i * (2 * pi / cfg.constants.lambda) * ...
        cfg.ris.dx * elementIndexX * path.xi(ell));
    steeringY = exp(1i * (2 * pi / cfg.constants.lambda) * ...
        cfg.ris.dy * elementIndexY * path.zeta(ell));
    steering = kron(steeringX, steeringY);

    for symbolIndex = 1:numSymbols
        spatialResponse = ...
            risReflections(:, symbolIndex, ell).' * steering;
        Y(:, symbolIndex) = Y(:, symbolIndex) + ...
            double(cfg.paths.visible(ell)) * ...
            path.rho(ell) * delayResponse * spatialResponse;
    end
end
end
