function [YSeparated, baseReflections, truth, generation] = ...
    separate_noiseless_observation(cfg)
%SEPARATE_NOISELESS_OBSERVATION Apply the historical DFT separation exactly.

[Y, truth, risReflections, codebook] = generate_noiseless_signal(cfg);
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

signalBlocks = reshape(Y, numSubcarriers, baseConfigCount, numBlocks);
YSeparated = complex(zeros(numSubcarriers, baseConfigCount, numRis));
for ell = 1:numRis
    for block = 1:numBlocks
        YSeparated(:, :, ell) = YSeparated(:, :, ell) + ...
            codebook(block, ell) * signalBlocks(:, :, block) / numBlocks;
    end
end

generation.rawNoiselessObservation = Y;
generation.risReflections = risReflections;
generation.codebook = codebook;
generation.separatedNoiseVariance = 0;
end
