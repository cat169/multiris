function result = b2_s4_original_sr( ...
    YSeparated, separatedNoiseVariance, baseReflections, ...
    positionPriorMessage, activeRisIndices, cfg)
%INFER_BASELINE2_FRAME_SYMMETRY_RESOLVED Fixed-posterior identity resolution.
% Channel estimation and Fisher weights are shared by at most two MAP calls.
activeRisIndices = activeRisIndices(:);
tChannel=tic; channel = estimate_channel_parameters( ...
    YSeparated, separatedNoiseVariance, baseReflections, activeRisIndices, cfg);
channelSeconds=toc(tChannel); tFim=tic; channelFimResults = cell(numel(activeRisIndices), 1);
for activeIndex = 1:numel(activeRisIndices)
    risIndex = activeRisIndices(activeIndex);
    qHat = [channel.xiEstimate(activeIndex); ...
        channel.zetaEstimate(activeIndex); ...
        cfg.ofdm.subcarrier_spacing * channel.delayEstimate(activeIndex)];
    channelFimResults{activeIndex} = compute_channel_fim( ...
        qHat, channel.rhoEstimate(activeIndex), baseReflections(:, :, risIndex), ...
        separatedNoiseVariance, cfg);
end
fimSeconds=toc(tFim); tMap=tic; [mapA, originalException] = b2_s4_try_map(channel, channelFimResults, ...
    positionPriorMessage, [], activeRisIndices, cfg);
mapSeconds=toc(tMap); solveCount=1; selected = mapA;
swapException=struct('identifier',"",'message',"");
diagnostic = struct('swapTriggered',false,'swapAccepted',false, ...
    'originalObjective',mapA.objective,'swappedInitialPriorCost',nan, ...
    'originalPriorCost',nan,'comparisonTolerance',nan,'refinedObjective',nan, ...
    'refinedExitflag',nan,'refinedIterations',nan, ...
    'originalPosition',mapA.positionEstimate,'originalExitflag',mapA.exitflag, ...
    'originalIterations',mapA.iterations,'swappedInitialPosition',nan(6,1), ...
    'refinedPosition',nan(6,1));
if mapA.exitflag > 0 && isfinite(mapA.objective) && all(isfinite(mapA.positionEstimate))
    xA = mapA.positionEstimate;
    xB0 = [xA(4:6); xA(1:3)];
    L = chol(positionPriorMessage.covariance, 'lower');
    m = positionPriorMessage.mean(:);
    rA = L \ (xA - m); rB = L \ (xB0 - m);
    priorCostA = 0.5 * (rA.' * rA); priorCostB = 0.5 * (rB.' * rB);
    tolerance = 100 * eps(max([1, abs(priorCostA), abs(priorCostB)]));
    diagnostic.originalPriorCost = priorCostA;
    diagnostic.swappedInitialPriorCost = priorCostB;
    diagnostic.comparisonTolerance = tolerance;
    diagnostic.swappedInitialPosition = xB0;
    if priorCostB < priorCostA - tolerance
        diagnostic.swapTriggered = true;
        tMap=tic; [mapB, swapException] = b2_s4_try_map(channel, channelFimResults, ...
            positionPriorMessage, xB0, activeRisIndices, cfg);
        mapSeconds=mapSeconds+toc(tMap); solveCount=solveCount+1; diagnostic.refinedObjective = mapB.objective;
        diagnostic.refinedExitflag = mapB.exitflag;
        diagnostic.refinedIterations = mapB.iterations;
        diagnostic.refinedPosition = mapB.positionEstimate;
        if mapB.exitflag > 0 && isfinite(mapB.objective) && ...
                all(isfinite(mapB.positionEstimate)) && mapB.objective < mapA.objective
            selected = mapB;
            diagnostic.swapAccepted = true;
        end
    end
end
result.originalException=originalException; result.swapException=swapException;
result.fisher=channelFimResults; result.channelSeconds=channelSeconds; result.fimSeconds=fimSeconds; result.mapSeconds=mapSeconds; result.solveCount=solveCount; result.channel = channel;
result.map = selected;
result.staticPositionEstimate = selected.staticPositionEstimate;
result.mobilePositionEstimate = selected.mobilePositionEstimate;
result.symmetryResolution = diagnostic;
end
