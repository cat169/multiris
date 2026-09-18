function result = run_paper_baseline_A( ...
    YSeparated, separatedNoiseVariance, baseReflections, cfg)
%RUN_PAPER_BASELINE_A Execute the authorized comparison-A algorithm.

activeRisIndices = cfg.paperBaselineA.oracleSupport(:);
searchInfo.lowerBound = cfg.paperBaselineA.lowerBound;
searchInfo.upperBound = cfg.paperBaselineA.upperBound;
searchInfo.gridStep = cfg.paperBaselineA.gridStep;

channel = estimate_channel_parameters(YSeparated, ...
    separatedNoiseVariance, baseReflections, activeRisIndices, cfg);
coarse = estimate_position_coarse(channel.delayEstimate, ...
    channel.xiEstimate, channel.zetaEstimate, activeRisIndices, ...
    searchInfo, cfg);
refined = refine_position_mle_eq57_identity(channel.delayEstimate, ...
    channel.xiEstimate, channel.zetaEstimate, activeRisIndices, ...
    coarse.staticPositionEstimate, coarse.mobilePositionEstimate, cfg);

result.activeRisIndices = activeRisIndices;
result.searchInfo = searchInfo;
result.channel = channel;
result.coarse = coarse;
result.refined = refined;
result.staticPositionEstimate = refined.staticPositionEstimate;
result.mobilePositionEstimate = refined.mobilePositionEstimate;
result.label = cfg.paperBaselineA.label;
end
