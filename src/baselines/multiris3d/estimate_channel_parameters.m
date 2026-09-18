function result = estimate_channel_parameters( ...
    YSeparated, separatedNoiseVariance, ...
    baseReflections, activeRisIndices, cfg)
%ESTIMATE_CHANNEL_PARAMETERS 集成 Multi-RIS 3D 信道参数估计流水线。
% activeRisIndices 是外部显式传入的 oracle support；噪声方差仅原样保留。

assert(size(YSeparated, 1) == cfg.ofdm.num_subcarriers);
assert(size(YSeparated, 3) == cfg.ris.count);
assert(size(baseReflections, 1) == cfg.ris.num_elements);
assert(size(baseReflections, 3) == cfg.ris.count);
assert(size(YSeparated, 2) == size(baseReflections, 2));
assert(~isempty(activeRisIndices));
assert(all(ismember(activeRisIndices, 1:cfg.ris.count)));
assert(isfinite(separatedNoiseVariance));
assert(separatedNoiseVariance >= 0);

activeRisIndices = activeRisIndices(:);

delayResult = estimate_delay_coarse( ...
    YSeparated, activeRisIndices, cfg);
spatialResult = estimate_spatial_coarse( ...
    YSeparated, baseReflections, ...
    activeRisIndices, delayResult.delayEstimate, cfg);
refinedResult = refine_channel_mle( ...
    YSeparated, baseReflections, activeRisIndices, ...
    delayResult.delayEstimate, spatialResult.xiEstimate, ...
    spatialResult.zetaEstimate, cfg);

result.activeRisIndices = activeRisIndices;
result.separatedNoiseVariance = separatedNoiseVariance;
result.delayEstimate = refinedResult.delayEstimate;
result.xiEstimate = refinedResult.xiEstimate;
result.zetaEstimate = refinedResult.zetaEstimate;
result.rhoEstimate = refinedResult.rhoEstimate;
result.coarseDelayEstimate = delayResult.delayEstimate;
result.coarseXiEstimate = spatialResult.xiEstimate;
result.coarseZetaEstimate = spatialResult.zetaEstimate;
result.initialResidualCost = refinedResult.initialResidualCost;
result.finalResidualCost = refinedResult.finalResidualCost;
result.initialProjectionScore = refinedResult.initialProjectionScore;
result.finalProjectionScore = refinedResult.finalProjectionScore;
result.exitflag = refinedResult.exitflag;
result.iterations = refinedResult.iterations;
result.funcCount = refinedResult.funcCount;
end
