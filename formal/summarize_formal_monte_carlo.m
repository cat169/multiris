function summary = summarize_formal_monte_carlo()
%SUMMARIZE_FORMAL_MONTE_CARLO Read checkpoints and write aggregate results.
% This function never calls signal generation or inference functions.

formalDirectory = fileparts(mfilename('fullpath'));
root = fileparts(formalDirectory);
checkpointDirectory = fullfile(root, 'formal_results', 'checkpoints');
outputDirectory = fullfile(root, 'formal_results', 'summary');
snrDbValues = (-20:5:30).';
expectedSeeds = (20261001:20261300).';
numSnr = numel(snrDbValues);

metricNames = {'coarseUE1Error','coarseUE2Error','coarseJointError', ...
    'refinedUE1Error','refinedUE2Error','refinedJointError'};
rmse = nan(numSnr, numel(metricNames));
medianError = nan(size(rmse));
p95Error = nan(size(rmse));
maximumError = nan(size(rmse));
metricSampleCount = zeros(size(rmse));
plannedCount = repmat(numel(expectedSeeds), numSnr, 1);
completedCount = zeros(numSnr, 1);
noOutputCount = zeros(numSnr, 1);
finiteNotConvergedCount = zeros(numSnr, 1);
normalReturnCount = zeros(numSnr, 1);
noOutputRate = nan(numSnr, 1);
finiteNotConvergedRate = nan(numSnr, 1);
finiteNotConvergedConditionalRate = nan(numSnr, 1);
runtimeMean = nan(numSnr, 1);
runtimeMedian = nan(numSnr, 1);
runtimeP95 = nan(numSnr, 1);
runtimeMaximum = nan(numSnr, 1);
runtimeSampleCount = zeros(numSnr, 1);
referenceFingerprint = [];

for snrIndex = 1:numSnr
    snrDb = snrDbValues(snrIndex);
    checkpointFile = fullfile(checkpointDirectory, ...
        ['checkpoint_snr_', snr_tag(snrDb), '.mat']);
    if ~isfile(checkpointFile)
        continue
    end
    loaded = load(checkpointFile, 'checkpoint');
    assert(isfield(loaded, 'checkpoint'), ...
        'paperA:InvalidCheckpoint', 'Missing checkpoint variable: %s', ...
        checkpointFile);
    checkpoint = loaded.checkpoint;
    validate_checkpoint_for_summary( ...
        checkpoint, snrDb, expectedSeeds, checkpointFile);
    if isempty(referenceFingerprint)
        referenceFingerprint = checkpoint.fingerprint;
    else
        assert(isequaln(referenceFingerprint, checkpoint.fingerprint), ...
            'paperA:MixedCheckpointVersions', ...
            'Checkpoint fingerprint differs: %s', checkpointFile);
    end

    completed = logical(checkpoint.completed(:));
    trials = checkpoint.trials(completed);
    completedCount(snrIndex) = numel(trials);
    if isempty(trials)
        continue
    end
    status = string({trials.status}).';
    noOutputCount(snrIndex) = sum(status == "1_no_valid_output");
    finiteNotConvergedCount(snrIndex) = sum( ...
        status == "2_finite_optimizer_not_converged");
    normalReturnCount(snrIndex) = sum(status == "3_finite_normal_return");
    noOutputRate(snrIndex) = noOutputCount(snrIndex) / numel(trials);
    finiteNotConvergedRate(snrIndex) = ...
        finiteNotConvergedCount(snrIndex) / numel(trials);
    finiteCount = sum([trials.hasRefinedFiniteOutput]);
    if finiteCount > 0
        finiteNotConvergedConditionalRate(snrIndex) = ...
            finiteNotConvergedCount(snrIndex) / finiteCount;
    end

    runtime = [trials.elapsedSeconds].';
    runtime = runtime(isfinite(runtime));
    runtimeSampleCount(snrIndex) = numel(runtime);
    if ~isempty(runtime)
        runtimeMean(snrIndex) = mean(runtime);
        runtimeMedian(snrIndex) = median(runtime);
        runtimeP95(snrIndex) = empirical_percentile(runtime, 95);
        runtimeMaximum(snrIndex) = max(runtime);
    end

    for metricIndex = 1:numel(metricNames)
        values = [trials.(metricNames{metricIndex})].';
        values = values(isfinite(values));
        metricSampleCount(snrIndex, metricIndex) = numel(values);
        if isempty(values)
            continue
        end
        rmse(snrIndex, metricIndex) = sqrt(mean(values.^2));
        medianError(snrIndex, metricIndex) = median(values);
        p95Error(snrIndex, metricIndex) = empirical_percentile(values, 95);
        maximumError(snrIndex, metricIndex) = max(values);
    end
end

assert(~isempty(referenceFingerprint), 'paperA:NoCheckpoints', ...
    'No formal checkpoints were found in %s.', checkpointDirectory);
summary.schemaVersion = 1;
summary.createdAt = char(datetime('now', ...
    'Format', 'yyyy-MM-dd HH:mm:ss Z'));
summary.conditioningNote = [ ...
    "Coarse statistics condition on finite coarse output. " + ...
    "Refined statistics condition on finite refined output; " + ...
    "finite nonconverged outputs are retained."];
summary.snrDbValues = snrDbValues;
summary.metricNames = metricNames;
summary.rmse = rmse;
summary.medianError = medianError;
summary.p95Error = p95Error;
summary.maximumError = maximumError;
summary.metricSampleCount = metricSampleCount;
summary.plannedCount = plannedCount;
summary.completedCount = completedCount;
summary.noOutputCount = noOutputCount;
summary.finiteNotConvergedCount = finiteNotConvergedCount;
summary.normalReturnCount = normalReturnCount;
summary.noOutputRate = noOutputRate;
summary.finiteNotConvergedRate = finiteNotConvergedRate;
summary.finiteNotConvergedConditionalRate = ...
    finiteNotConvergedConditionalRate;
summary.runtimeMean = runtimeMean;
summary.runtimeMedian = runtimeMedian;
summary.runtimeP95 = runtimeP95;
summary.runtimeMaximum = runtimeMaximum;
summary.runtimeSampleCount = runtimeSampleCount;
summary.fingerprint = referenceFingerprint;

summaryTable = build_summary_table(summary);
if ~isfolder(outputDirectory)
    mkdir(outputDirectory);
end
save(fullfile(outputDirectory, 'formal_summary.mat'), ...
    'summary', 'summaryTable', '-v7.3');
writetable(summaryTable, ...
    fullfile(outputDirectory, 'formal_summary.csv'));
end

function tableValue = build_summary_table(summary)
tableValue = table(summary.snrDbValues, summary.plannedCount, ...
    summary.completedCount, summary.noOutputCount, ...
    summary.noOutputRate, summary.finiteNotConvergedCount, ...
    summary.finiteNotConvergedRate, ...
    summary.finiteNotConvergedConditionalRate, ...
    summary.normalReturnCount, summary.runtimeSampleCount, ...
    summary.runtimeMean, summary.runtimeMedian, summary.runtimeP95, ...
    summary.runtimeMaximum, 'VariableNames', { ...
    'snrDb','plannedCount','completedCount','noOutputCount', ...
    'noOutputRate','finiteNotConvergedCount', ...
    'finiteNotConvergedRate','finiteNotConvergedConditionalRate', ...
    'normalReturnCount','runtimeSampleCount','runtimeMean', ...
    'runtimeMedian','runtimeP95','runtimeMaximum'});
for metricIndex = 1:numel(summary.metricNames)
    stem = summary.metricNames{metricIndex};
    tableValue.([stem, '_rmse']) = summary.rmse(:, metricIndex);
    tableValue.([stem, '_median']) = ...
        summary.medianError(:, metricIndex);
    tableValue.([stem, '_p95']) = summary.p95Error(:, metricIndex);
    tableValue.([stem, '_maximum']) = ...
        summary.maximumError(:, metricIndex);
    tableValue.([stem, '_sampleCount']) = ...
        summary.metricSampleCount(:, metricIndex);
end
end

function validate_checkpoint_for_summary( ...
    checkpoint, snrDb, expectedSeeds, checkpointFile)
required = {'schemaVersion','fingerprint','snrDb','seeds', ...
    'numTrials','completed','trials'};
assert(all(isfield(checkpoint, required)), ...
    'paperA:CheckpointFieldsMismatch', ...
    'Incomplete checkpoint: %s', checkpointFile);
assert(checkpoint.schemaVersion == 1 && checkpoint.snrDb == snrDb, ...
    'paperA:CheckpointConditionMismatch', ...
    'Schema or SNR mismatch: %s', checkpointFile);
assert(checkpoint.numTrials == numel(expectedSeeds) && ...
    isequal(checkpoint.seeds(:), expectedSeeds), ...
    'paperA:CheckpointSeedMismatch', ...
    'Trial count or seed mismatch: %s', checkpointFile);
assert(numel(checkpoint.completed) == numel(expectedSeeds) && ...
    numel(checkpoint.trials) == numel(expectedSeeds), ...
    'paperA:CheckpointSizeMismatch', ...
    'Checkpoint array sizes mismatch: %s', checkpointFile);
completedIndices = find(checkpoint.completed(:));
for index = completedIndices.'
    trial = checkpoint.trials(index);
    assert(trial.trialIndex == index && ...
        trial.seed == expectedSeeds(index) && trial.snrDb == snrDb, ...
        'paperA:CheckpointTrialIdentityMismatch', ...
        'Completed trial identity mismatch in %s at index %d.', ...
        checkpointFile, index);
end
end

function value = empirical_percentile(values, percentile)
values = sort(values(:));
index = max(1, ceil(percentile / 100 * numel(values)));
value = values(index);
end

function tag = snr_tag(snrDb)
if snrDb < 0
    tag = sprintf('m%03d', abs(snrDb));
else
    tag = sprintf('p%03d', snrDb);
end
end
