function run_baseline2_s4_equivalence(trackingRoot, outputDirectory)
% Read-only 1600-result audit, followed by 13 original-seed paired recomputes.
% Does not run a new MC stream, modify old results, or invoke VB/B1/PEB.
here=fileparts(mfilename('fullpath')); repo=fileparts(fileparts(here));
if nargin<2,outputDirectory=fullfile(repo,'test_results/baseline2_s4_interface');end
if ~isfolder(outputDirectory),mkdir(outputDirectory);end
oldPath=path;oldRng=rng;cleanup=onCleanup(@()restore(oldPath,oldRng));
inference=fullfile(repo,'baseline2/inference');
standalone={inference,fullfile(repo,'src/baselines/multiris3d'),fullfile(repo,'src/structured/current/signal')};
current=fullfile(trackingRoot,'src/structured/current');
history=fullfile(trackingRoot,'reports/baseline2_multistart_independent_validation');
C=load(fullfile(history,'checkpoint.mat'),'records','completed','spec','files','signature');
assert(numel(C.records)==1600&&all(C.completed));
% Fixed before any recomputation; selection tolerance remains the original 1e-10.
tolerances=struct('positionInf',1e-8,'objectiveRelative',1e-10, ...
    'initialInf',0,'fisherMaxAbs',0,'gradientRelative',1e-8);
audit=cell(1600,1);groups=strings(1600,1);
for i=1:1600
    z=C.records{i};q=z.candidates(1:4);L=chol(z.prior.covariance,'lower');
    initial=z.prior.mean;
    if z.srDiagnostics.swapAccepted,initial=z.srDiagnostics.swappedInitialPosition;end
    X=[initial z.prior.mean z.prior.mean+L(:,1) z.prior.mean-L(:,1)];
    assert(isequal(string({q.label}),["originalSR","mean","+TxX","-TxX"]));
    assert(isequal([q.initial],X));
    mask=[q.finite]&[q.objectiveConsistent];indices=find(mask);index=0;
    if ~isempty(indices),[~,k]=min([q(indices).objective]);index=indices(k);end
    assert(index==z.selectedA(2)&&index>0);
    assert(isequal(q(index).position,z.candidates(z.selectedA(2)).position));
    finiteCheck=arrayfun(@(a)isreal(a.position)&&all(isfinite(a.position))&& ...
        isreal(a.reportedObjective)&&isfinite(a.reportedObjective),q);
    assert(isequal(finiteCheck(:).',[q.finite]));
    assert(all(strlength(string({q.exceptionIdentifier}))==0));
    positive=find(mask&([q.exitflag]>0));indexB=0;
    if ~isempty(positive),[~,k]=min([q(positive).objective]);indexB=positive(k);end
    assert(indexB==z.selectedB(2));
    groups(i)=z.job.pattern+"/"+string(z.job.snrDb);
    selected=q(index).position;
    row=struct('jobIndex',i,'pattern',z.job.pattern,'snrDb',z.job.snrDb,'trial',z.job.trial, ...
        'selectedIndex',index,'positiveExitflagIndex',indexB,'finite',true, ...
        'expectedMapCalls',z.originalSolves+3,'objective',q(index).objective, ...
        'txX',selected(1),'txY',selected(2),'txZ',selected(3), ...
        'rxX',selected(4),'rxY',selected(5),'rxZ',selected(6));
    for j=1:4
        row.(sprintf('candidate%dObjective',j))=q(j).objective;
        row.(sprintf('candidate%dExitflag',j))=q(j).exitflag;
    end
    audit{i}=row;
end
audit=struct2table([audit{:}]);writetable(audit,fullfile(outputDirectory,'retained_1600_audit.csv'));
fprintf('RETAINED_AUDIT_PASSED 1600/1600; this is a saved-result audit, not rerunning 1600.\n');
% One trial from every originally validated pattern/SNR group; both repaired
% 011/10 dB cases and the -15 dB SR-swap/accuracy-degradation case are mandatory.
[~,first]=unique(groups,'stable');jobs=first(:).';
for i=1:1600
    j=C.records{i}.job;
    if (j.pattern=="011"&&j.snrDb==10&&ismember(j.trial,[27 209])) || ...
            C.records{i}.srDiagnostics.swapTriggered
        jobs(end+1)=i; %#ok<AGROW>
    end
end
jobs=unique(jobs,'stable');assert(numel(jobs)==13);
rows=cell(numel(jobs),1);pairedDetails=cell(numel(jobs),1);
for n=1:numel(jobs)
    z=C.records{jobs(n)};job=z.job;cfg=C.spec.configuration;
    cfg.paths.visible=(char(job.pattern).'-'0')==1;cfg.noise.snr_db=job.snrDb;
    restoredefaultpath;
    addpath(current,fullfile(current,'signal'), ...
        fullfile(current,'experiment/fixed_noise_single_frame'), ...
        fullfile(repo,'src/baselines/multiris3d'),inference, ...
        fullfile(trackingRoot,'analysis/baseline2_multistart_validation'), ...
        fullfile(trackingRoot,'analysis/baseline2_multistart_independent_validation'));
    rng(job.signalSeed,'twister');
    [Y,nv,W]=generate_fixed_noise_observation(cfg,C.spec.amplitudes(job.snrIndex),C.spec.rawNoiseVariance);
    assert(isequal(nv,z.noiseVariance));
    original=evaluate_fixed_candidates(Y,nv,W,z.prior,z.active,cfg);
    assert(isequaln(original.channel,z.channel),'baseline2S4:ReplayMismatch','Original MLE differs from retained record.');
    assert(isequaln(original.fisher,z.fisher));
    % The new call has no Tracking path, including no temporary analysis path.
    restoredefaultpath;addpath(standalone{:});
    assert(isempty(which('evaluate_fixed_candidates'))&&isempty(which('cached_sr')));
    assert(strcmpi(which('upa_steering'),fullfile(repo,'src/structured/current/signal/upa_steering.m')));
    profile clear;profile on;
    actual=infer_baseline2_frame_s4(Y,nv,W,z.prior,z.active,cfg);
    profile off;p=profile('info');
    counts=[profileCalls(p,'estimate_channel_parameters'),profileCalls(p,'compute_channel_fim'),profileCalls(p,'refine_position_map')];
    expected=[1 numel(z.active) z.originalSolves+3];assert(isequal(counts,expected));
    assert(isequaln(actual.channel,original.channel)&&isequaln(actual.fisher,original.fisher));
    maxPosition=0;maxObjective=0;maxGradient=0;
    for j=1:4
        a=actual.candidates(j);b=original.candidates(j);saved=z.candidates(j);
        assert(isequal(a.initial,b.initial)&&isequal(a.initial,saved.initial)&&a.label==b.label);
        assert(a.finite==b.finite&&a.objectiveConsistent==b.objectiveConsistent);
        assert(a.exitflag==b.exitflag&&a.iterations==b.iterations);
        assert(a.exceptionIdentifier==b.exceptionIdentifier&&a.exceptionMessage==b.exceptionMessage);
        dx=max(norm(a.position-b.position,inf),norm(a.position-saved.position,inf));
        dj=max(abs(a.objective-b.objective),abs(a.objective-saved.objective));
        dg=abs(a.gradientInf-b.gradientInf);
        assert(dx<=tolerances.positionInf&&dj<=tolerances.objectiveRelative*max(1,abs(b.objective)));
        assert(dg<=tolerances.gradientRelative*max(1,abs(b.gradientInf)));
        maxPosition=max(maxPosition,dx);maxObjective=max(maxObjective,dj);maxGradient=max(maxGradient,dg);
    end
    assert(actual.selectedIndex==original.selectedA(2)&&actual.selectedIndex==z.selectedA(2));
    assert(actual.positiveExitflagIndex==original.selectedB(2));
    assert(isequal(actual.positionEstimate,original.candidates(original.selectedA(2)).position));
    rows{n}=struct('jobIndex',jobs(n),'pattern',job.pattern,'snrDb',job.snrDb,'trial',job.trial, ...
        'maxPositionDifference',maxPosition,'maxObjectiveDifference',maxObjective, ...
        'maxGradientDifference',maxGradient,'selectedIndex',actual.selectedIndex, ...
        'channelMleCalls',counts(1),'fisherCalls',counts(2),'mapCalls',counts(3), ...
        'swapTriggered',actual.symmetryResolution.swapTriggered,'finite',actual.finite, ...
        'jointError',norm(actual.positionEstimate-[cfg.geometry.static_user;cfg.geometry.mobile_user]));
    pairedDetails{n}=struct('job',job,'actual',actual,'referenceFour',original.candidates(1:4));
    fprintf('S4_RECOMPUTE %d/%d %s %+g trial=%d dx=%.12g dJ=%.12g calls=%d\n', ...
        n,numel(jobs),job.pattern,job.snrDb,job.trial,maxPosition,maxObjective,counts(3));
end
recompute=struct2table([rows{:}]);writetable(recompute,fullfile(outputDirectory,'representative_recompute.csv'));
save(fullfile(outputDirectory,'equivalence_summary.mat'),'audit','recompute','pairedDetails','tolerances','jobs','-v7');
fprintf('S4_EQUIVALENCE_PASSED retained=1600 recomputed=%d maxDx=%.15g maxDJ=%.15g\n', ...
    height(recompute),max(recompute.maxPositionDifference),max(recompute.maxObjectiveDifference));
end
function count=profileCalls(p,name)
items=p.FunctionTable;mask=strcmp({items.FunctionName},name);
assert(nnz(mask)==1,'baseline2S4:ProfileName','Expected profiler entry %s.',name);
count=items(mask).NumCalls;
end
function restore(p,r)
profile off;path(p);rng(r);
end
