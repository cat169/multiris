function result=b1_pai_evaluate_candidates(initializations,channel,cfg)
% Selection is solely finite real position + original finite objective.
template=struct('initialization',struct(),'estimate',nan(6,1),'initialObjective',nan,'finalObjective',nan,'spatialObjective',nan,'delayObjective',nan,'exitflag',nan,'iterations',nan,'funcCount',nan,'firstOrderOptimality',nan,'minimumRisDistance',nan,'stepFromInitial',nan,'elapsedSeconds',nan,'valid',false,'qualityFailed',true,'strictGradient',false,'exceptionIdentifier','','exceptionMessage','','formal',struct());
candidates=repmat(template,1,numel(initializations)); geometry=cfg.geometry; geometry=rmfield(geometry,{'static_user','mobile_user'});
for k=1:numel(initializations)
    q=template; init=initializations(k); q.initialization=init; q.initialObjective=init.initialObjective; timer=tic;
    try
        f=refine_position_mle_eq57_identity(channel.delayEstimate,channel.xiEstimate,channel.zetaEstimate,channel.activeRisIndices,init.x0(1:3),init.x0(4:6),cfg);
    catch ex
        % Only known candidate-level numerical failures are recoverable.
        numericalIds={'refine_position_mle_eq57_identity:PositionAtRisCenter', ...
            'optimlib:commonMsgs:UsrObjUndefAtX0','optimlib:commonMsgs:UserJacUndefAtX0'};
        if ~any(strcmp(ex.identifier,numericalIds)), rethrow(ex); end
        q.exceptionIdentifier=ex.identifier; q.exceptionMessage=ex.message;
        q.elapsedSeconds=toc(timer); candidates(k)=q; continue;
    end
        q.formal=f; q.estimate=[f.staticPositionEstimate;f.mobilePositionEstimate]; [q.finalObjective,parts,q.minimumRisDistance]=b1_pai_objective(q.estimate,channel,geometry,cfg.constants.c);
        q.spatialObjective=parts(1); q.delayObjective=parts(2); q.exitflag=f.exitflag; q.iterations=f.iterations; q.funcCount=f.funcCount; q.firstOrderOptimality=f.firstOrderOptimality; q.stepFromInitial=norm(q.estimate-init.x0);
        validateattributes(q.estimate,{'numeric'},{'numel',6});
        validateattributes(f.initialCost,{'numeric'},{'scalar'});
        validateattributes(f.finalCost,{'numeric'},{'scalar'});
        validateattributes(q.initialObjective,{'numeric'},{'scalar'});
        validateattributes(q.finalObjective,{'numeric'},{'scalar'});
        finiteValues=[isreal(q.estimate)&&all(isfinite(q.estimate)), ...
            isreal(f.initialCost)&&isfinite(f.initialCost),isreal(f.finalCost)&&isfinite(f.finalCost), ...
            isreal(q.initialObjective)&&isfinite(q.initialObjective),isreal(q.finalObjective)&&isfinite(q.finalObjective)];
        q.valid=all(finiteValues);
        if ~q.valid
            labels={'position','optimizer initial objective','optimizer final objective','independent initial objective','independent final objective'};
            q.exceptionIdentifier='baseline1PAI:NonfiniteCandidate';
            q.exceptionMessage=['Nonfinite or nonreal candidate values: ',strjoin(labels(~finiteValues),', '),'.'];
            q.elapsedSeconds=toc(timer); candidates(k)=q; continue;
        end
        assert(abs(q.initialObjective-f.initialCost)<=1e-12*max(1,f.initialCost),'fallback:InitialObjectiveMismatch','Initial objective differs from formal residual.');
        assert(abs(q.finalObjective-f.finalCost)<=1e-12*max(1,f.finalCost),'fallback:FinalObjectiveMismatch','Final objective differs from formal residual.');
        q.qualityFailed=~q.valid || f.exitflag<=0 || ~isfinite(f.firstOrderOptimality) || f.firstOrderOptimality>1e-6 || q.minimumRisDistance<=1e-6 || q.finalObjective>q.initialObjective+1e-10*max(1,q.initialObjective);
        q.strictGradient=q.valid && f.firstOrderOptimality<=1e-10;
    q.elapsedSeconds=toc(timer); candidates(k)=q;
end
valid=find([candidates.valid]); selectedIndex=0;
if ~isempty(valid), [~,j]=min([candidates(valid).finalObjective]); selectedIndex=valid(j); end
result=struct('candidates',candidates,'selectedIndex',selectedIndex,'finite',selectedIndex>0,'estimate',nan(6,1),'solveCount',numel(candidates),'elapsedSolverSeconds',sum([candidates.elapsedSeconds]));
if selectedIndex>0, result.estimate=candidates(selectedIndex).estimate; end
end
