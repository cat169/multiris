function result = estimate_delay_coarse( ...
    YSeparated, activeRisIndices, cfg)
%ESTIMATE_DELAY_COARSE 估计 active RIS 路径的粗时延。
% YSeparated 为 numSubcarriers × baseConfigCount × numRis 的逐 RIS
% 分离观测。activeRisIndices 是外部显式传入的 oracle path support；本函数
% 不估计路径存在变量，也不读取 cfg.paths.visible 来决定 support。

assert(size(YSeparated, 1) == cfg.ofdm.num_subcarriers);
assert(size(YSeparated, 3) == cfg.ris.count);
assert(~isempty(activeRisIndices));
assert(all(ismember(activeRisIndices, 1:cfg.ris.count)));

% 论文复现采用固定的 1024 点 coarse delay grid。当前 steering convention
% 的共轭为正指数，因此沿第一维使用 IFFT。
nFft = 1024;
subcarrierSpacing = cfg.ofdm.subcarrier_spacing;
delayGrid = (0:(nFft - 1)).' / (nFft * subcarrierSpacing);

activeRisIndices = activeRisIndices(:);
numActiveRis = numel(activeRisIndices);
delayEstimate = zeros(numActiveRis, 1);
delayBin = zeros(numActiveRis, 1);
delayScore = zeros(nFft, numActiveRis);

for activeIndex = 1:numActiveRis
    ell = activeRisIndices(activeIndex);
    Hell = YSeparated(:, :, ell);
    delayTransform = ifft(Hell, nFft, 1);
    delayScore(:, activeIndex) = sum(abs(delayTransform).^2, 2);

    [~, peakIndex] = max(delayScore(:, activeIndex));
    delayBin(activeIndex) = peakIndex - 1;
    delayEstimate(activeIndex) = ...
        delayBin(activeIndex) / (nFft * subcarrierSpacing);
end

result.activeRisIndices = activeRisIndices;
result.delayEstimate = delayEstimate;
result.delayGrid = delayGrid;
result.delayScore = delayScore;
result.delayBin = delayBin;
result.nFft = nFft;
end
