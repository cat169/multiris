function run_b1_pai_hardening_validation(trackingRoot,outputDirectory)
% Fixed frozen-data regression; no new MLE, VB, Baseline2 or Monte Carlo.
here=fileparts(mfilename('fullpath')); repo=fileparts(fileparts(here));
if nargin<2,outputDirectory=fullfile(repo,'test_results','baseline1_pai_hardening');end
if ~isfolder(outputDirectory),mkdir(outputDirectory);end
oldPath=path;oldRng=rng;cleanup=onCleanup(@() restore(oldPath,oldRng));
addpath(fullfile(trackingRoot,'src/structured/current'),'-begin');addpath(fullfile(trackingRoot,'src/structured/current/signal'),'-begin');addpath(fullfile(trackingRoot,'src/structured/current/experiment/fixed_noise_single_frame'),'-begin');addpath(fullfile(repo,'src/baselines/multiris3d'),'-begin');addpath(here,'-begin');
clear infer_baseline1_frame_projection_assisted estimate_position_coarse refine_position_mle_eq57_identity estimate_channel_parameters;
paths={'infer_baseline1_frame_projection_assisted','estimate_channel_parameters','estimate_position_coarse','refine_position_mle_eq57_identity'};
fid=fopen(fullfile(outputDirectory,'matlab_path_check.txt'),'w');fprintf(fid,'MATLAB %s\n',version);
for k=1:numel(paths),resolved=which(paths{k});assert(startsWith(strrep(resolved,'\','/'),strrep(repo,'\','/')));fprintf(fid,'%s: %s\n',paths{k},resolved);end
fclose(fid);
production={which('infer_baseline1_frame_projection_assisted'),fullfile(repo,'src/baselines/multiris3d/private/b1_pai_evaluate_candidates.m'),fullfile(repo,'src/baselines/multiris3d/private/b1_pai_build_initializations.m')};
[dependencies,products]=matlab.codetools.requiredFilesAndProducts(production);save(fullfile(outputDirectory,'dependency_check.mat'),'dependencies','products');
analysis=cell(size(production));for k=1:numel(production),analysis{k}=checkcode(production{k},'-id');end
save(fullfile(outputDirectory,'code_analysis.mat'),'production','analysis');
S=load(fullfile(trackingRoot,'reports/baseline1_pai_validation/checkpoint.mat'),'checkpoint');previous=S.checkpoint;
H=load(fullfile(trackingRoot,'reports/baseline1_fallback_comparison/checkpoint.mat'),'checkpoint');raw=H.checkpoint;
rows={};errors=[];savedCandidates={};
for i=1:numel(previous.historical)
    old=previous.historical{i};input=raw.records{i};cfg=raw.datasets(input.job.datasetIndex).configuration;
    new=infer_baseline1_frame_projection_assisted(input.YSeparated,input.separatedNoiseVariance,input.baseReflections,cfg,input.channel);
    stage='historical_normal';if old.job.failure,stage='historical_failure';end
    row=compare_pair(old.pai,new,stage,old.job.pattern,-15,old.job.trial,input.signalSeed);
    rows{end+1}=row;savedCandidates{end+1}=struct('stage',stage,'pattern',old.job.pattern,'trial',old.job.trial,'result',new); %#ok<AGROW>
    if old.job.failure,errors(end+1)=norm([new.staticPositionEstimate;new.mobilePositionEstimate]-[cfg.geometry.static_user;cfg.geometry.mobile_user]);end %#ok<AGROW>
    fprintf('REGRESSION %d/262 %s %s trial %d delta=%.12g\n',i,stage,old.job.pattern,old.job.trial,row.positionDifference);
end
for i=1:numel(previous.high)
    old=previous.high{i};d=previous.datasets(old.job.datasetIndex);t=d.highTrials(old.job.highTrialIndex);cfg=d.configuration;cfg.noise.snr_db=old.job.snrDb;
    rng(old.signalSeed,'twister');[Y,nv,W]=generate_fixed_noise_observation(cfg,d.amplitudes(old.job.snrIndex),d.settings.rawNoiseVariance);
    new=infer_baseline1_frame_projection_assisted(Y,nv,W,cfg,old.pai.channel);
    row=compare_pair(old.pai,new,'high_normal',old.job.pattern,old.job.snrDb,old.job.trial,old.signalSeed);
    assert(~new.fallbackTriggered&&new.extraPositionSolves==0);
    rows{end+1}=row;savedCandidates{end+1}=struct('stage','high_normal','pattern',old.job.pattern,'trial',old.job.trial,'result',new); %#ok<AGROW>
    fprintf('REGRESSION %d/262 high_normal %s %g trial %d delta=%.12g\n',202+i,old.job.pattern,old.job.snrDb,old.job.trial,row.positionDifference);
end
rmse=sqrt(mean(errors.^2));assert(numel(errors)==162&&abs(rmse-61.3577904760193)<=1e-6);
regression=struct2table([rows{:}]);writetable(regression,fullfile(outputDirectory,'regression.csv'));save(fullfile(outputDirectory,'regression_candidates.mat'),'savedCandidates','rmse');
results=runtests(fullfile(here,'Baseline1PaiHardeningTest.m'));disp(results);assert(all([results.Passed]));
boundary=table(string({results.Name}).',[results.Passed].',[results.Failed].',[results.Incomplete].',[results.Duration].','VariableNames',{'name','passed','failed','incomplete','durationSeconds'});writetable(boundary,fullfile(outputDirectory,'boundary_results.csv'));
fid=fopen(fullfile(outputDirectory,'validation_summary.txt'),'w');fprintf(fid,'Historical failures: 162/162 pass\nHistorical normals: 40/40 pass\nHigh-SNR normals: 60/60 pass\nHistorical Joint RMSE: %.15g m\nBoundary tests: %d/%d pass\nMLE reused for all regressions; no MLE inference or new Monte Carlo.\n',rmse,nnz([results.Passed]),numel(results));fclose(fid);
fprintf('HARDENING_VALIDATION_PASSED rmse=%.15g boundary=%d\n',rmse,numel(results));
end
function row=compare_pair(old,new,stage,pattern,snr,trial,seed)
assert(old.finite==new.finite&&old.fallbackTriggered==new.fallbackTriggered&&old.selectedIndex==new.selectedIndex&&old.candidateCount==new.candidateCount&&old.extraPositionSolves==new.extraPositionSolves&&isequaln(old.channel,new.channel));
delta=norm([old.staticPositionEstimate;old.mobilePositionEstimate]-[new.staticPositionEstimate;new.mobilePositionEstimate],inf);
assert(delta<=1e-8&&old.qualityFailed==new.qualityFailed&&isequaln(old.qualityReasons,new.qualityReasons));
maxObjective=0;maxCandidatePosition=0;flagDifferences=0;iterationDifferences=0;qualityDifferences=0;
for k=1:numel(old.candidates)
    a=old.candidates(k);b=new.candidates(k);position=norm(a.estimate-b.estimate,inf);objective=abs(a.finalObjective-b.finalObjective);
    assert(position<=1e-8&&objective<=1e-12*max(1,abs(a.finalObjective))&&abs(a.initialObjective-b.initialObjective)<=1e-12*max(1,abs(a.initialObjective)));
    assert(a.valid==b.valid&&a.exitflag==b.exitflag&&a.iterations==b.iterations&&a.funcCount==b.funcCount&&a.qualityFailed==b.qualityFailed&&a.strictGradient==b.strictGradient&&isequaln(a.qualityReasons,b.qualityReasons));
    assert(isequaln(a.initialization,b.initialization)&&abs(a.firstOrderOptimality-b.firstOrderOptimality)<=1e-8*max(1,abs(a.firstOrderOptimality))&&abs(a.minimumRisDistance-b.minimumRisDistance)<=1e-8);
    maxObjective=max(maxObjective,objective);maxCandidatePosition=max(maxCandidatePosition,position);flagDifferences=flagDifferences+(a.exitflag~=b.exitflag);iterationDifferences=iterationDifferences+(a.iterations~=b.iterations);qualityDifferences=qualityDifferences+(a.qualityFailed~=b.qualityFailed);
end
if ~new.fallbackTriggered,assert(isequaln(old.coarse,new.coarse)&&new.extraPositionSolves==0);end
row=struct('stage',stage,'pattern',pattern,'snrDb',snr,'trial',trial,'signalSeed',seed,'finite',new.finite,'fallbackTriggered',new.fallbackTriggered,'selectedIndex',new.selectedIndex,'candidateCount',new.candidateCount,'positionDifference',delta,'maxCandidatePositionDifference',maxCandidatePosition,'maxObjectiveDifference',maxObjective,'exitflagDifferences',flagDifferences,'iterationDifferences',iterationDifferences,'qualityDifferences',qualityDifferences,'mleEqual',isequaln(old.channel,new.channel),'extraPositionSolves',new.extraPositionSolves,'passed',true);
end
function restore(p,r),path(p);rng(r);end
