function result = estimate_spatial_coarse( ...
    YSeparated, baseReflections, ...
    activeRisIndices, delayEstimate, cfg)
%ESTIMATE_SPATIAL_COARSE 使用 coherent matched filter 粗估计空间频率。
% activeRisIndices 和对应的 delayEstimate 必须由外部显式传入；本函数不读取
% 路径存在真值来决定 support，也不使用真实 delay 或空间参数。

assert(size(YSeparated, 1) == cfg.ofdm.num_subcarriers);
assert(size(YSeparated, 3) == cfg.ris.count);
assert(size(baseReflections, 1) == cfg.ris.num_elements);
assert(size(baseReflections, 3) == cfg.ris.count);
assert(size(YSeparated, 2) == size(baseReflections, 2));
assert(~isempty(activeRisIndices));
assert(all(ismember(activeRisIndices, 1:cfg.ris.count)));
assert(numel(delayEstimate) == numel(activeRisIndices));

activeRisIndices = activeRisIndices(:);
delayEstimate = delayEstimate(:);
numActiveRis = numel(activeRisIndices);

gridStep = 0.02;
xiGrid = (-1:gridStep:(1 - gridStep)).';
zetaGrid = (-1:gridStep:(1 - gridStep)).';
subcarrierIndex = (0:(cfg.ofdm.num_subcarriers - 1)).';

xiEstimate = zeros(numActiveRis, 1);
zetaEstimate = zeros(numActiveRis, 1);
% 索引约定：dimension 1 -> xi，dimension 2 -> zeta，
% dimension 3 -> active RIS（而不是全体 RIS 的原始索引）。
spatialScore = zeros( ...
    numel(xiGrid), numel(zetaGrid), numActiveRis);

for activeIndex = 1:numActiveRis
    ell = activeRisIndices(activeIndex);
    Hell = YSeparated(:, :, ell);
    W = baseReflections(:, :, ell);
    tauHat = delayEstimate(activeIndex);
    delayVector = exp(-1i * 2 * pi * subcarrierIndex * ...
        cfg.ofdm.subcarrier_spacing * tauHat);

    for ix = 1:numel(xiGrid)
        for iz = 1:numel(zetaGrid)
            a = upa_steering(xiGrid(ix), zetaGrid(iz), cfg);
            s = W.' * a;
            U = delayVector * s.';
            spatialScore(ix, iz, activeIndex) = ...
                abs(sum(U .* conj(Hell), 'all'))^2;
        end
    end

    scoreCurrent = spatialScore(:, :, activeIndex);
    [~, linearIndex] = max(scoreCurrent(:));
    [ix, iz] = ind2sub(size(scoreCurrent), linearIndex);
    xiEstimate(activeIndex) = xiGrid(ix);
    zetaEstimate(activeIndex) = zetaGrid(iz);
end

assert(all(isfinite(spatialScore(:))));

result.activeRisIndices = activeRisIndices;
result.xiEstimate = xiEstimate;
result.zetaEstimate = zetaEstimate;
result.xiGrid = xiGrid;
result.zetaGrid = zetaGrid;
result.spatialScore = spatialScore;
result.gridStep = gridStep;
end
