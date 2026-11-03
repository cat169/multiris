function check_baseline2_s4_standalone(repo, inputFile)
% One frozen-input check with only this checkout's three production paths.
S=load(inputFile,'input');input=S.input;
oldPath=path;cleanup=onCleanup(@()path(oldPath));
restoredefaultpath;
addpath(fullfile(repo,'baseline2/inference'), ...
    fullfile(repo,'src/baselines/multiris3d'),fullfile(repo,'src/structured/current/signal'));
names={'infer_baseline2_frame_s4','estimate_channel_parameters','compute_channel_fim', ...
    'refine_position_map','upa_steering'};
for k=1:numel(names)
    assert(startsWith(which(names{k}),repo),'baseline2S4:ExternalDependency',names{k});
end
assert(isempty(which('cached_sr'))&&isempty(which('evaluate_fixed_candidates')));
[dependencies,products]=matlab.codetools.requiredFilesAndProducts(which(names{1}));
for k=1:numel(dependencies)
    assert(~contains(dependencies{k},'Tracking_code'),'baseline2S4:TrackingDependency',dependencies{k});
end
result=infer_baseline2_frame_s4(input.Y,input.nv,input.W,input.prior,input.active,input.cfg);
assert(isequal(result.positionEstimate,input.referencePosition));
assert(result.selectedIndex==input.referenceIndex&&result.callCounts.map==input.expectedMapCalls);
fprintf('STANDALONE_CHECK_PASSED checkout=%s dependencies=%d mapCalls=%d products=%d\n', ...
    repo,numel(dependencies),result.callCounts.map,numel(products));
end
