function result=infer_baseline1_frame_pai_c(YSeparated,separatedNoiseVariance,baseReflections,cfg,channel)
%INFER_BASELINE1_FRAME_PAI_C Independent distance-scaled B1-PAI-C candidate.
% Four arguments match run_paper_baseline_A. Optional fifth argument reuses
% an already computed ORIGINAL MLE on these same observations for pairing.
% Only estimate_position_coarse:NoValidCandidate triggers projection.
% All C objectives use xi/zeta and delay in metres. Original B1 is intact.
% The original MLE, coarse geometry, grid and ranker are reused unchanged.
timer=tic; channelTimer=tic;
if nargin<5 || isempty(channel)
    channel=estimate_channel_parameters(YSeparated,separatedNoiseVariance,baseReflections,cfg.paperBaselineA.oracleSupport(:),cfg);
    channelReused=false;
else
    validate_external_channel(channel,cfg.paperBaselineA.oracleSupport,separatedNoiseVariance);
    channelReused=true;
end
channelSeconds=toc(channelTimer); active=cfg.paperBaselineA.oracleSupport(:);
search=struct('lowerBound',cfg.paperBaselineA.lowerBound,'upperBound',cfg.paperBaselineA.upperBound,'gridStep',cfg.paperBaselineA.gridStep);
coarse=struct(); originalExceptionIdentifier=''; originalExceptionMessage=''; coarseTimer=tic;
try
    coarse=estimate_position_coarse(channel.delayEstimate,channel.xiEstimate,channel.zetaEstimate,active,search,cfg);
catch ex
    if ~strcmp(ex.identifier,'estimate_position_coarse:NoValidCandidate'), rethrow(ex); end
    originalExceptionIdentifier=ex.identifier; originalExceptionMessage=ex.message;
end
coarseSeconds=toc(coarseTimer); fallbackTriggered=~isempty(originalExceptionIdentifier); originalRefinementSeconds=0; fallbackSeconds=0; geometry=rmfield(cfg.geometry,{'static_user','mobile_user'});
if ~fallbackTriggered
    refinementTimer=tic;
    refined=refine_position_mle_eq57_distance_scaled(channel.delayEstimate,channel.xiEstimate,channel.zetaEstimate,active,coarse.staticPositionEstimate,coarse.mobilePositionEstimate,cfg);
    originalRefinementSeconds=toc(refinementTimer);
    x=[refined.staticPositionEstimate;refined.mobilePositionEstimate]; [cost,parts,minDistance]=b1_pai_objective_c(x,channel,geometry,cfg.constants.c);
    init=struct('x0',[coarse.staticPositionEstimate;coarse.mobilePositionEstimate],'gridIndex',coarse.bestGridIndex,'risIndex',0,'initialObjective',refined.initialCost);
    candidate=struct('initialization',init,'estimate',x,'initialObjective',refined.initialCost,'finalObjective',cost,'spatialObjective',parts(1),'delayObjective',parts(2),'exitflag',refined.exitflag,'iterations',refined.iterations,'funcCount',refined.funcCount,'firstOrderOptimality',refined.firstOrderOptimality,'minimumRisDistance',minDistance,'stepFromInitial',norm(x-init.x0),'elapsedSeconds',originalRefinementSeconds,'valid',isreal(x)&&all(isfinite(x))&&isreal(cost)&&isfinite(cost),'qualityFailed',false,'strictGradient',false,'exceptionIdentifier','','exceptionMessage','','formal',refined);
    candidates=candidate; selectedIndex=1; initializationSource='original_geometry'; initializationDiagnostics=struct(); extraSolves=0;
else
    fallbackTimer=tic;
    [initializations,initializationDiagnostics]=b1_pai_build_initializations_c(channel,geometry,search,cfg.constants.c);
    evaluation=b1_pai_evaluate_candidates_c(initializations,channel,cfg);
    candidates=evaluation.candidates; selectedIndex=evaluation.selectedIndex; extraSolves=evaluation.solveCount; fallbackSeconds=toc(fallbackTimer); initializationSource='projection_assisted';
end
for k=1:numel(candidates)
    candidates(k).qualityReasons=quality_reasons(candidates(k));
    candidates(k).qualityFailed=any(structfun(@(v) logical(v),candidates(k).qualityReasons));
    candidates(k).strictGradient=candidates(k).valid && isfinite(candidates(k).firstOrderOptimality) && candidates(k).firstOrderOptimality<=1e-10;
end
result=struct('activeRisIndices',active,'searchInfo',search,'channel',channel,'coarse',coarse,'refined',struct(),'staticPositionEstimate',nan(3,1),'mobilePositionEstimate',nan(3,1),'label','Baseline 1-PAI-C (Distance-Scaled Delay Residual)','fallbackTriggered',fallbackTriggered,'initializationSource',initializationSource,'initializationDiagnostics',initializationDiagnostics,'candidates',candidates,'candidateCount',numel(candidates),'selectedIndex',selectedIndex,'finite',false,'qualityFailed',true,'qualityReasons',struct(),'strictGradient',false,'initialObjective',nan,'finalObjective',nan,'originalExceptionIdentifier',originalExceptionIdentifier,'originalExceptionMessage',originalExceptionMessage,'exceptionIdentifier','','exceptionMessage','','channelReused',channelReused,'extraPositionSolves',extraSolves);
if selectedIndex>0
    chosen=candidates(selectedIndex); result.staticPositionEstimate=chosen.estimate(1:3); result.mobilePositionEstimate=chosen.estimate(4:6); result.refined=chosen.formal;
    result.finite=chosen.valid; result.qualityFailed=chosen.qualityFailed; result.qualityReasons=chosen.qualityReasons; result.strictGradient=chosen.strictGradient; result.initialObjective=chosen.initialObjective; result.finalObjective=chosen.finalObjective;
else
    result.exceptionIdentifier='baseline1PAI:NoValidRefinement'; result.exceptionMessage='No finite position with a computable distance-scaled objective from projected initializations.';
end
result.objectiveDefinition='sum(deltaXi.^2+deltaZeta.^2+distanceResidual.^2)';
result.delayObjectiveUnits='metres_squared';
result.timing=struct('channelSeconds',channelSeconds,'originalCoarseSeconds',coarseSeconds,'originalRefinementSeconds',originalRefinementSeconds,'fallbackSeconds',fallbackSeconds,'totalSeconds',toc(timer));
end
function reasons=quality_reasons(q)
reasons=struct('nonpositiveExitflag',~isfinite(q.exitflag)||q.exitflag<=0,'nonfiniteFirstOrder',~isfinite(q.firstOrderOptimality),'largeFirstOrder',isfinite(q.firstOrderOptimality)&&q.firstOrderOptimality>1e-6,'nearRisCenter',isfinite(q.minimumRisDistance)&&q.minimumRisDistance<=1e-6,'objectiveIncrease',isfinite(q.finalObjective)&&q.finalObjective>q.initialObjective+1e-10*max(1,q.initialObjective),'nonfinitePosition',~isreal(q.estimate)||any(~isfinite(q.estimate)),'nonfiniteObjective',~isreal(q.finalObjective)||~isfinite(q.finalObjective));
end
function validate_external_channel(channel,support,noiseVariance)
% Structural consistency only; the caller guarantees observation provenance.
required={'activeRisIndices','delayEstimate','xiEstimate','zetaEstimate'};
assert(isstruct(channel)&&isscalar(channel)&&all(isfield(channel,required)), ...
    'baseline1PAI:InvalidExternalChannel','External channel must be a scalar structure with the required MLE fields.');
active=channel.activeRisIndices;
assert(isnumeric(active)&&isvector(active)&&isreal(active)&&all(isfinite(active(:)))&&isequal(active(:),support(:)), ...
    'baseline1PAI:ExternalSupportMismatch','External MLE support must match oracleSupport, including order.');
for k=2:numel(required)
    value=channel.(required{k});
    assert(isnumeric(value)&&isvector(value)&&numel(value)==numel(support), ...
        'baseline1PAI:ExternalChannelSize','Each external MLE estimate must be a vector with one entry per active RIS.');
    assert(isreal(value)&&all(isfinite(value(:))), ...
        'baseline1PAI:ExternalChannelNonfinite','External delay and spatial frequency estimates must be finite real values.');
end
if isfield(channel,'separatedNoiseVariance')
    value=channel.separatedNoiseVariance;
    assert(isnumeric(value)&&isscalar(value)&&isreal(value)&&isfinite(value)&&isequal(value,noiseVariance), ...
        'baseline1PAI:ExternalNoiseMismatch','External MLE noise variance must match the current input.');
end
end
