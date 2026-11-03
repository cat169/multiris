function audit_frozen_independent_guards(trackingRoot,outputDirectory)
% Read saved 800-case evidence only. Never invokes MLE or position solvers.
S=load(fullfile(trackingRoot,'reports/baseline1_pai_validation/checkpoint.mat'),'checkpoint');cp=S.checkpoint;rows={};
allowed={'refine_position_mle_eq57_identity:PositionAtRisCenter','optimlib:commonMsgs:UsrObjUndefAtX0','optimlib:commonMsgs:UserJacUndefAtX0'};
for i=1:numel(cp.independent)
    r=cp.independent{i};cfg=cp.datasets(r.job.datasetIndex).configuration;p=r.pai;ch=p.channel;support=cfg.paperBaselineA.oracleSupport;
    channelValid=isstruct(ch)&&isscalar(ch)&&isequal(ch.activeRisIndices(:),support(:));
    for name={'delayEstimate','xiEstimate','zetaEstimate'},v=ch.(name{1});channelValid=channelValid&&isnumeric(v)&&isvector(v)&&numel(v)==numel(support)&&isreal(v)&&all(isfinite(v));end
    noiseValid=~isfield(ch,'separatedNoiseVariance')||isequal(ch.separatedNoiseVariance,r.separatedNoiseVariance);
    a=cfg.paperBaselineA;[x,y,z]=ndgrid(a.lowerBound(1):a.gridStep:a.upperBound(1),a.lowerBound(2):a.gridStep:a.upperBound(2),a.lowerBound(3):a.gridStep:a.upperBound(3));grid=[x(:).';y(:).';z(:).'];minimumDistance=inf;
    for l=support(:).',minimumDistance=min(minimumDistance,min(vecnorm(grid-cfg.geometry.ris_centers(:,l))));end
    nonfiniteCandidates=0;newFatalExceptions=0;
    for k=1:numel(p.candidates)
        q=p.candidates(k);
        if ~isempty(q.exceptionIdentifier),newFatalExceptions=newFatalExceptions+~any(strcmp(q.exceptionIdentifier,allowed));continue;end
        f=q.formal;finite=isreal(q.estimate)&&all(isfinite(q.estimate))&&isreal(q.initialObjective)&&isfinite(q.initialObjective)&&isreal(q.finalObjective)&&isfinite(q.finalObjective)&&isreal(f.initialCost)&&isfinite(f.initialCost)&&isreal(f.finalCost)&&isfinite(f.finalCost);
        nonfiniteCandidates=nonfiniteCandidates+~finite;
    end
    row=struct('pattern',r.job.pattern,'snrDb',r.job.snrDb,'trial',r.job.trial,'externalChannelValid',channelValid,'noiseMatches',noiseValid,'minimumTxGridRisDistance',minimumDistance,'candidateCount',p.candidateCount,'nonfiniteCandidates',nonfiniteCandidates,'previouslySwallowedNowFatal',newFatalExceptions);
    assert(channelValid&&noiseValid&&minimumDistance>0&&nonfiniteCandidates==0&&newFatalExceptions==0,'paiHardening:FrozenGuardChanged','Saved independent evidence enters a changed branch.');rows{end+1}=row; %#ok<AGROW>
end
writetable(struct2table([rows{:}]),fullfile(outputDirectory,'frozen_800_guard_audit.csv'));
fprintf('FROZEN_GUARD_AUDIT_PASSED trials=%d candidates=%d; no optimization\n',numel(rows),sum(cellfun(@(r) r.candidateCount,rows)));
end
