function run_formal_monte_carlo()
%RUN_FORMAL_MONTE_CARLO Run or resume the fixed comparison-A experiment.
% This entry performs 11 SNR conditions x 300 independent single-frame
% trials. It never reads or writes the validation results directory.

originalPath = path;
originalRng = rng;
stateCleanup = onCleanup(@() restore_state(originalPath, originalRng)); %#ok<NASGU>
formalDirectory = fileparts(mfilename('fullpath'));
root = fileparts(formalDirectory);
restoredefaultpath;
addpath(fullfile(root, 'src', 'structured', 'current'), '-begin');
addpath(fullfile(root, 'src', 'structured', 'current', 'signal'), '-begin');
addpath(fullfile(root, 'src', 'baselines', 'multiris3d'), '-begin');
addpath(fullfile(root, 'experiment'), '-begin');
addpath(formalDirectory, '-begin');

cfg = paper_baseline_A_config();
snrDbValues = (-20:5:30).';
seeds = (20261001:20261300).';
assert(numel(snrDbValues) == 11);
assert(numel(seeds) == 300);
assert(isempty(intersect(seeds, (20260918:20260922).')));

fingerprint = build_fingerprint(root, cfg, snrDbValues, seeds);
resultsDirectory = fullfile(root, 'formal_results');
checkpointDirectory = fullfile(resultsDirectory, 'checkpoints');
if ~isfolder(checkpointDirectory)
    mkdir(checkpointDirectory);
end

for snrIndex = 1:numel(snrDbValues)
    snrDb = snrDbValues(snrIndex);
    conditionCfg = cfg;
    conditionCfg.noise.snr_db = snrDb;
    checkpointFile = fullfile(checkpointDirectory, ...
        ['checkpoint_snr_', snr_tag(snrDb), '.mat']);
    checkpoint = load_or_initialize_checkpoint( ...
        checkpointFile, fingerprint, conditionCfg, snrDb, seeds);
    fprintf('\nFormal comparison A: SNR %g dB, %d/%d already complete.\n', ...
        snrDb, nnz(checkpoint.completed), numel(seeds));

    for trialIndex = 1:numel(seeds)
        if checkpoint.completed(trialIndex)
            continue
        end
        trial = run_one_trial( ...
            conditionCfg, snrDb, seeds(trialIndex), trialIndex);
        checkpoint.trials(trialIndex) = trial;
        checkpoint.completed(trialIndex) = true;
        checkpoint.updatedAt = char(datetime('now', ...
            'Format', 'yyyy-MM-dd HH:mm:ss Z'));
        atomic_save_checkpoint(checkpointFile, checkpoint);
        fprintf(['SNR=%g trial=%d/%d seed=%d status=%s ', ...
            'coarseJoint=%.6g refinedJoint=%.6g time=%.3fs\n'], ...
            snrDb, trialIndex, numel(seeds), trial.seed, trial.status, ...
            trial.coarseJointError, trial.refinedJointError, ...
            trial.elapsedSeconds);
    end
end
fprintf('\nAll formal comparison-A checkpoints are complete.\n');
end

function trial = run_one_trial(cfgTemplate, snrDb, seed, trialIndex)
cfg = cfgTemplate;
cfg.noise.snr_db = snrDb;
rng(seed, 'twister');
trial = empty_trial();
trial.trialIndex = trialIndex;
trial.seed = seed;
trial.snrDb = snrDb;
timer = tic;

try
    trial.failureStage = "signal_generation";
    [YSeparated, separatedNoiseVariance, baseReflections] = ...
        generate_noisy_separated_signal(cfg);
    trial.separatedNoiseVariance = separatedNoiseVariance;

    trial.failureStage = "channel_estimation";
    active = cfg.paperBaselineA.oracleSupport(:);
    channel = estimate_channel_parameters(YSeparated, ...
        separatedNoiseVariance, baseReflections, active, cfg);
    trial.channel = channel;
    trial.hasChannelOutput = true;

    trial.failureStage = "coarse_position";
    searchInfo.lowerBound = cfg.paperBaselineA.lowerBound;
    searchInfo.upperBound = cfg.paperBaselineA.upperBound;
    searchInfo.gridStep = cfg.paperBaselineA.gridStep;
    coarse = estimate_position_coarse(channel.delayEstimate, ...
        channel.xiEstimate, channel.zetaEstimate, active, ...
        searchInfo, cfg);
    trial.coarse = coarse;
    trial.hasCoarseOutput = all(isfinite([ ...
        coarse.staticPositionEstimate; coarse.mobilePositionEstimate]));
    if trial.hasCoarseOutput
        trial.coarseUE1Position = coarse.staticPositionEstimate;
        trial.coarseUE2Position = coarse.mobilePositionEstimate;
        trial.coarseUE1ErrorVector = ...
            coarse.staticPositionEstimate - cfg.geometry.static_user;
        trial.coarseUE2ErrorVector = ...
            coarse.mobilePositionEstimate - cfg.geometry.mobile_user;
        [trial.coarseUE1Error, trial.coarseUE2Error, ...
            trial.coarseJointError] = position_errors( ...
            coarse.staticPositionEstimate, coarse.mobilePositionEstimate, cfg);
        trial.numValidCandidates = coarse.numValidCandidates;
    end

    trial.failureStage = "position_refinement";
    refined = refine_position_mle_eq57_identity( ...
        channel.delayEstimate, channel.xiEstimate, ...
        channel.zetaEstimate, active, ...
        coarse.staticPositionEstimate, coarse.mobilePositionEstimate, cfg);
    trial.refined = refined;
    trial.hasRefinedFiniteOutput = all(isfinite([ ...
        refined.staticPositionEstimate; refined.mobilePositionEstimate]));
    if trial.hasRefinedFiniteOutput
        trial.refinedUE1Position = refined.staticPositionEstimate;
        trial.refinedUE2Position = refined.mobilePositionEstimate;
        trial.refinedUE1ErrorVector = ...
            refined.staticPositionEstimate - cfg.geometry.static_user;
        trial.refinedUE2ErrorVector = ...
            refined.mobilePositionEstimate - cfg.geometry.mobile_user;
        [trial.refinedUE1Error, trial.refinedUE2Error, ...
            trial.refinedJointError] = position_errors( ...
            refined.staticPositionEstimate, refined.mobilePositionEstimate, cfg);
        trial.positionUpdateNorm = norm([ ...
            refined.staticPositionEstimate - coarse.staticPositionEstimate; ...
            refined.mobilePositionEstimate - coarse.mobilePositionEstimate]);
    end
    trial.positionInitialObjective = refined.initialCost;
    trial.positionFinalObjective = refined.finalCost;
    trial.positionExitflag = refined.exitflag;
    trial.positionIterations = refined.iterations;
    trial.positionFirstOrderOptimality = refined.firstOrderOptimality;
    trial.channelOptimizerConverged = all(channel.exitflag > 0);
    trial.positionOptimizerConverged = refined.exitflag > 0;

    if ~trial.hasRefinedFiniteOutput
        trial.status = "1_no_valid_output";
    elseif ~(trial.channelOptimizerConverged && ...
            trial.positionOptimizerConverged)
        trial.status = "2_finite_optimizer_not_converged";
    else
        trial.status = "3_finite_normal_return";
    end
    trial.failureStage = "";
catch exception
    trial.exceptionIdentifier = string(exception.identifier);
    trial.exceptionMessage = string(exception.message);
    trial.status = "1_no_valid_output";
end
trial.elapsedSeconds = toc(timer);
end

function checkpoint = load_or_initialize_checkpoint( ...
    checkpointFile, fingerprint, cfg, snrDb, seeds)
if isfile(checkpointFile)
    loaded = load(checkpointFile, 'checkpoint');
    assert(isfield(loaded, 'checkpoint') && isscalar(loaded.checkpoint), ...
        'paperA:InvalidCheckpoint', ...
        'Checkpoint does not contain one scalar checkpoint: %s', ...
        checkpointFile);
    checkpoint = loaded.checkpoint;
    validate_checkpoint(checkpoint, fingerprint, cfg, snrDb, seeds, ...
        checkpointFile);
    return
end

checkpoint.schemaVersion = 1;
checkpoint.experimentLabel = cfg.paperBaselineA.label;
checkpoint.createdAt = char(datetime('now', ...
    'Format', 'yyyy-MM-dd HH:mm:ss Z'));
checkpoint.updatedAt = checkpoint.createdAt;
checkpoint.fingerprint = fingerprint;
checkpoint.cfg = cfg;
checkpoint.snrDb = snrDb;
checkpoint.seeds = seeds;
checkpoint.numTrials = numel(seeds);
checkpoint.completed = false(numel(seeds), 1);
checkpoint.trials = repmat(empty_trial(), numel(seeds), 1);
end

function validate_checkpoint( ...
    checkpoint, fingerprint, cfg, snrDb, seeds, checkpointFile)
required = {'schemaVersion', 'fingerprint', 'cfg', 'snrDb', 'seeds', ...
    'numTrials', 'completed', 'trials'};
assert(all(isfield(checkpoint, required)), ...
    'paperA:CheckpointFieldsMismatch', ...
    'Checkpoint fields are incomplete: %s', checkpointFile);
assert(checkpoint.schemaVersion == 1, ...
    'paperA:CheckpointSchemaMismatch', ...
    'Checkpoint schema mismatch: %s', checkpointFile);
assert(isequaln(checkpoint.fingerprint, fingerprint), ...
    'paperA:CheckpointFingerprintMismatch', ...
    'Configuration or source hashes differ from checkpoint: %s', ...
    checkpointFile);
assert(isequaln(checkpoint.cfg, cfg), ...
    'paperA:CheckpointConfigMismatch', ...
    'Full configuration differs from checkpoint: %s', checkpointFile);
assert(isequal(checkpoint.snrDb, snrDb), ...
    'paperA:CheckpointSnrMismatch', ...
    'SNR differs from checkpoint: %s', checkpointFile);
assert(isequal(checkpoint.seeds(:), seeds(:)) && ...
    checkpoint.numTrials == numel(seeds), ...
    'paperA:CheckpointSeedMismatch', ...
    'Trial count or seeds differ from checkpoint: %s', checkpointFile);
assert(numel(checkpoint.completed) == numel(seeds) && ...
    numel(checkpoint.trials) == numel(seeds), ...
    'paperA:CheckpointSizeMismatch', ...
    'Checkpoint array sizes differ: %s', checkpointFile);
end

function fingerprint = build_fingerprint(root, cfg, snrDbValues, seeds)
relativeFiles = { ...
    fullfile('src','structured','current','simulation_config.m'); ...
    fullfile('src','structured','current','signal','generate_noiseless_signal.m'); ...
    fullfile('src','structured','current','signal','generate_noisy_separated_signal.m'); ...
    fullfile('src','structured','current','signal','upa_steering.m'); ...
    fullfile('src','baselines','multiris3d','estimate_delay_coarse.m'); ...
    fullfile('src','baselines','multiris3d','estimate_spatial_coarse.m'); ...
    fullfile('src','baselines','multiris3d','refine_channel_mle.m'); ...
    fullfile('src','baselines','multiris3d','estimate_channel_parameters.m'); ...
    fullfile('src','baselines','multiris3d','estimate_position_coarse.m'); ...
    fullfile('src','baselines','multiris3d','refine_position_mle_eq57_identity.m'); ...
    fullfile('experiment','paper_baseline_A_config.m'); ...
    fullfile('experiment','run_paper_baseline_A.m'); ...
    fullfile('formal','run_formal_monte_carlo.m')};
hashes = strings(numel(relativeFiles), 1);
for index = 1:numel(relativeFiles)
    absoluteFile = fullfile(root, relativeFiles{index});
    assert(isfile(absoluteFile), 'paperA:MissingSourceFile', ...
        'Required source file is missing: %s', absoluteFile);
    hashes(index) = sha256_file(absoluteFile);
end
fingerprint.schemaVersion = 1;
fingerprint.snrDbValues = snrDbValues;
fingerprint.seeds = seeds;
fingerprint.numTrialsPerSnr = numel(seeds);
fingerprint.sourceFiles = relativeFiles;
fingerprint.sourceSha256 = hashes;
fingerprint.visibility = cfg.paths.visible;
fingerprint.oracleSupport = cfg.paperBaselineA.oracleSupport;
fingerprint.searchLowerBound = cfg.paperBaselineA.lowerBound;
fingerprint.searchUpperBound = cfg.paperBaselineA.upperBound;
fingerprint.searchGridStep = cfg.paperBaselineA.gridStep;
end

function hash = sha256_file(fileName)
digest = javaMethod('getInstance', ...
    'java.security.MessageDigest', 'SHA-256');
fid = fopen(fileName, 'rb');
assert(fid >= 0, 'paperA:HashOpenFailed', ...
    'Cannot open source file for hashing: %s', fileName);
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
while ~feof(fid)
    bytes = fread(fid, 1024 * 1024, '*uint8');
    if ~isempty(bytes)
        digest.update(bytes);
    end
end
bytes = typecast(digest.digest(), 'uint8');
hash = lower(string(reshape(dec2hex(bytes, 2).', 1, [])));
end

function atomic_save_checkpoint(checkpointFile, checkpoint)
checkpointDirectory = fileparts(checkpointFile);
temporaryFile = [checkpointFile, '.tmp'];
if isfile(temporaryFile)
    delete(temporaryFile);
end
try
    save(temporaryFile, 'checkpoint', '-v7.3');
    [moved, message] = movefile(temporaryFile, checkpointFile, 'f');
    if ~moved
        error('paperA:CheckpointMoveFailed', ...
            'Cannot replace checkpoint %s: %s', checkpointFile, message);
    end
catch exception
    if isfile(temporaryFile)
        delete(temporaryFile);
    end
    error('paperA:CheckpointSaveFailed', ...
        'Checkpoint save failed in %s: %s', ...
        checkpointDirectory, exception.message);
end
end

function [ue1Error, ue2Error, jointError] = position_errors(p1, p2, cfg)
ue1Error = norm(p1 - cfg.geometry.static_user);
ue2Error = norm(p2 - cfg.geometry.mobile_user);
jointError = hypot(ue1Error, ue2Error);
end

function tag = snr_tag(snrDb)
if snrDb < 0
    tag = sprintf('m%03d', abs(snrDb));
else
    tag = sprintf('p%03d', snrDb);
end
end

function trial = empty_trial()
trial = struct( ...
    'trialIndex', nan, 'seed', nan, 'snrDb', nan, ...
    'status', "", 'failureStage', "", ...
    'exceptionIdentifier', "", 'exceptionMessage', "", ...
    'elapsedSeconds', nan, 'separatedNoiseVariance', nan, ...
    'hasChannelOutput', false, 'hasCoarseOutput', false, ...
    'hasRefinedFiniteOutput', false, ...
    'channelOptimizerConverged', false, ...
    'positionOptimizerConverged', false, ...
    'channel', struct(), 'coarse', struct(), 'refined', struct(), ...
    'coarseUE1Position', nan(3,1), 'coarseUE2Position', nan(3,1), ...
    'coarseUE1ErrorVector', nan(3,1), ...
    'coarseUE2ErrorVector', nan(3,1), ...
    'refinedUE1Position', nan(3,1), 'refinedUE2Position', nan(3,1), ...
    'refinedUE1ErrorVector', nan(3,1), ...
    'refinedUE2ErrorVector', nan(3,1), ...
    'coarseUE1Error', nan, 'coarseUE2Error', nan, ...
    'coarseJointError', nan, 'refinedUE1Error', nan, ...
    'refinedUE2Error', nan, 'refinedJointError', nan, ...
    'numValidCandidates', nan, 'positionUpdateNorm', nan, ...
    'positionInitialObjective', nan, 'positionFinalObjective', nan, ...
    'positionExitflag', nan, 'positionIterations', nan, ...
    'positionFirstOrderOptimality', nan);
end

function restore_state(originalPath, originalRng)
path(originalPath);
rng(originalRng);
end
