function result = infer_baseline2_frame_s4(Y, nv, W, prior, active, cfg)
%INFER_BASELINE2_FRAME_S4 Original SR plus the three validated fixed starts.
% Same six inputs as B2-SR. Selection is the validated finite-consistent A
% rule; positive exitflag is reported separately and never filters that rule.
% Required paths: baseline2/inference, src/baselines/multiris3d and
% src/structured/current/signal, all from this repository.
timer = tic;
active = active(:);
assert(isequal(size(prior.mean),[6,1]) && ...
    isequal(size(prior.covariance),[6,6]), ...
    'baseline2S4:InvalidPrior','Use a six-dimensional column mean and covariance.');
sr = b2_s4_original_sr(Y,nv,W,prior,active,cfg);
% Preserve the analysis evaluator's post-SR factorization and evaluation order.
L = chol(prior.covariance,'lower');
roots = cell(numel(active),1);
for j=1:numel(active)
    F=sr.fisher{j}.fim; F=(F+F.')/2; [V,D]=eig(F); e=real(diag(D));
    assert(all(e>=-100*eps*max(norm(F,2),realmin)));
    e(e<0)=0; roots{j}=diag(sqrt(e))*V.';
end
X=[prior.mean prior.mean+L(:,1) prior.mean-L(:,1)];
labels=["mean","+TxX","-TxX"];
candidates=repmat(candidateTemplate(),4,1);
x0=prior.mean;
if sr.symmetryResolution.swapAccepted
    x0=sr.symmetryResolution.swappedInitialPosition;
end
candidates(1)=assess(sr.map,x0,"originalSR",sr.mapSeconds, ...
    sr.originalException.identifier,sr.originalException.message);
for j=1:3
    t=tic;
    [map,exception]=b2_s4_try_map(sr.channel,sr.fisher,prior,X(:,j),active,cfg);
    candidates(j+1)=assess(map,X(:,j),labels(j),toc(t),exception.identifier,exception.message);
end
eligible=[candidates.finite]&[candidates.objectiveConsistent];
selectedIndex=pick(eligible);
positiveExitflagIndex=pick(eligible&([candidates.exitflag]>0));
result=struct('channel',sr.channel,'fisher',{sr.fisher}, ...
    'symmetryResolution',sr.symmetryResolution,'candidates',candidates, ...
    'selectedIndex',selectedIndex,'positiveExitflagIndex',positiveExitflagIndex, ...
    'selectionRule',"A: finite real position/objective and objective-consistent minimum full J", ...
    'finite',selectedIndex>0,'converged',false,'candidateStatus',"", ...
    'exceptionIdentifier',"",'exceptionMessage',"", ...
    'staticPositionEstimate',nan(3,1),'mobilePositionEstimate',nan(3,1), ...
    'positionEstimate',nan(6,1),'map',struct(),'objective',nan, ...
    'srOriginalException',sr.originalException,'srSwapException',sr.swapException, ...
    'callCounts',struct('channelMle',1,'fisher',numel(active),'map',sr.solveCount+3), ...
    'timing',struct('channelSeconds',sr.channelSeconds,'fisherSeconds',sr.fimSeconds, ...
    'mapSeconds',sum([candidates.seconds]),'totalSeconds',toc(timer)));
if all(eligible)
    result.candidateStatus="allCandidatesValid";
elseif any(eligible)
    result.candidateStatus="partialCandidateFailure";
else
    result.candidateStatus="allCandidatesFailed";
    result.exceptionIdentifier="baseline2S4:NoValidCandidate";
    result.exceptionMessage="No finite, real and objective-consistent candidate among the validated four.";
end
if selectedIndex>0
    chosen=candidates(selectedIndex);
    result.positionEstimate=chosen.position;
    result.staticPositionEstimate=chosen.position(1:3);
    result.mobilePositionEstimate=chosen.position(4:6);
    result.map=chosen.map;
    result.objective=chosen.objective;
    result.converged=chosen.exitflag>0;
end

    function c=assess(map,initial,label,seconds,id,msg)
        c=candidateTemplate(); c.label=label; c.initial=initial;
        c.position=map.positionEstimate; c.reportedObjective=map.objective;
        c.exitflag=map.exitflag; c.iterations=map.iterations; c.seconds=seconds;
        c.exceptionIdentifier=string(id); c.exceptionMessage=string(msg); c.map=map;
        c.finite=isreal(c.position)&&all(isfinite(c.position))&& ...
            isreal(map.objective)&&isfinite(map.objective);
        if c.finite
            try
                [res,jac]=b2_s4_residual(c.position,sr.channel,roots,prior.mean,L\eye(6),active,cfg);
            catch exception
                if ~strcmp(exception.identifier,'refine_position_map:PositionAtRisCenter')
                    rethrow(exception);
                end
                c.exceptionIdentifier=string(exception.identifier);
                c.exceptionMessage=string(exception.message);
                return;
            end
            n=3*numel(active);
            c.fisherObjective=.5*sum(res(1:n).^2);
            c.priorObjective=.5*sum(res(n+(1:6)).^2);
            c.objective=.5*(res.'*res); c.gradientInf=norm(jac.'*res,inf);
            c.objectiveConsistent=abs(c.objective-map.objective)<=1e-10*max(1,abs(map.objective));
        end
    end
    function idx=pick(mask)
        indices=find(mask); idx=0;
        if ~isempty(indices)
            [~,bestOffset]=min([candidates(indices).objective]); idx=indices(bestOffset);
        end
    end
end
function c=candidateTemplate()
c=struct('label',"",'initial',nan(6,1),'position',nan(6,1), ...
    'reportedObjective',nan,'fisherObjective',nan,'priorObjective',nan, ...
    'objective',nan,'gradientInf',nan,'exitflag',nan,'iterations',nan, ...
    'seconds',nan,'finite',false,'objectiveConsistent',false, ...
    'exceptionIdentifier',"",'exceptionMessage',"",'map',struct());
end
