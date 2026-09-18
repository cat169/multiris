function validation = run_all_validations()
%RUN_ALL_VALIDATIONS Run the authorized A/B/C validation sequence only.

originalPath = path;
originalRng = rng;
stateCleanup = onCleanup(@() restore_state(originalPath, originalRng)); %#ok<NASGU>
root = fileparts(fileparts(mfilename('fullpath')));
restoredefaultpath;
addpath(fullfile(root, 'src', 'structured', 'current'), '-begin');
addpath(fullfile(root, 'src', 'structured', 'current', 'signal'), '-begin');
addpath(fullfile(root, 'src', 'baselines', 'multiris3d'), '-begin');
addpath(fullfile(root, 'experiment'), '-begin');
addpath(fullfile(root, 'validation'), '-begin');

validation.createdAt = char(datetime('now', ...
    'Format', 'yyyy-MM-dd HH:mm:ss Z'));
validation.matlabVersion = version;
validation.label = [ ...
    "Multi-RIS 3D adapted to the single-frame two-user scenario, " + ...
    "with oracle support and a truth-centered coarse search region."];
validation.static = static_check();

cfg = paper_baseline_A_config();
resultsRoot = fullfile(root, 'results');
if ~isfolder(resultsRoot), mkdir(resultsRoot); end

validation.A = validate_truth_geometry(cfg);
save(fullfile(resultsRoot, 'validation_A_truth_geometry.mat'), ...
    'validation', '-v7.3');
assert(validation.A.structuralPass, 'paperA:StructuralFailureA', ...
    'Validation A found a structural error; B/C were not run.');

validation.B = validate_noiseless_pipeline(cfg);
save(fullfile(resultsRoot, 'validation_AB_noiseless.mat'), ...
    'validation', '-v7.3');
assert(validation.B.structuralPass, 'paperA:StructuralFailureB', ...
    'Validation B found a structural error; C was not run.');

validation.C = validate_small_noisy(cfg, resultsRoot);
save(fullfile(resultsRoot, 'validation_all.mat'), 'validation', '-v7.3');
writetable(struct2table(validation.C.summary), ...
    fullfile(resultsRoot, 'validation_C_summary.csv'));
write_text_summary(validation, fullfile(resultsRoot, 'validation_report.txt'));
end

function result = validate_truth_geometry(cfg)
[~, truth] = generate_noiseless_signal(cfg);
active = cfg.paperBaselineA.oracleSupport;
search = search_info(cfg);
coarse = estimate_position_coarse(truth.tau(active), truth.xi(active), ...
    truth.zeta(active), active, search, cfg);
refined = refine_position_mle_eq57_identity(truth.tau(active), ...
    truth.xi(active), truth.zeta(active), active, ...
    coarse.staticPositionEstimate, coarse.mobilePositionEstimate, cfg);

result.truth = truth;
result.principalXi = periodic_wrap(truth.xi(active), 2);
result.principalZeta = periodic_wrap(truth.zeta(active), 2);
result.normalProjectionUE1 = normal_projections( ...
    cfg.geometry.static_user, active, cfg);
result.normalProjectionUE2 = normal_projections( ...
    cfg.geometry.mobile_user, active, cfg);
result.searchInfo = search;
result.expectedGridSize = [11; 11; 11];
result.expectedCandidateCount = 1331;
result.coarse = coarse;
result.refined = refined;
result.coarseErrors = position_errors(coarse.staticPositionEstimate, ...
    coarse.mobilePositionEstimate, cfg);
result.refinedErrors = position_errors(refined.staticPositionEstimate, ...
    refined.mobilePositionEstimate, cfg);
result.residualComponents = residual_components(refined.finalResidual);
result.jacobianCheck = finite_difference_jacobian_check( ...
    coarse.staticPositionEstimate, coarse.mobilePositionEstimate, ...
    truth.tau(active), truth.xi(active), truth.zeta(active), ...
    active, cfg, refined.initialJacobian);
result.structuralPass = coarse.numGridCandidates == 1331 && ...
    isequal(coarse.gridSize, [11; 11; 11]) && ...
    coarse.numValidCandidates > 0 && ...
    all(isfinite([refined.staticPositionEstimate; ...
    refined.mobilePositionEstimate])) && ...
    result.jacobianCheck.angleMaxAbsError < 1e-6 && ...
    result.jacobianCheck.delayMaxAbsError < 1e-12;
end

function result = validate_noiseless_pipeline(cfg)
[YSeparated, baseReflections, truth, generation] = ...
    separate_noiseless_observation(cfg);
pipeline = run_paper_baseline_A(YSeparated, 0, baseReflections, cfg);
active = pipeline.activeRisIndices;
reconstruction = reconstruct_separated(pipeline.channel, ...
    baseReflections, active, cfg);

result.truth = truth;
result.generation = generation;
result.YSeparated = YSeparated;
result.baseReflections = baseReflections;
result.pipeline = pipeline;
result.channelErrors = channel_errors(pipeline.channel, truth, active);
difference = YSeparated(:, :, active) - reconstruction(:, :, active);
reference = YSeparated(:, :, active);
result.reconstructionRelativeResidual = norm(difference(:)) / ...
    norm(reference(:));
result.coarseErrors = position_errors( ...
    pipeline.coarse.staticPositionEstimate, ...
    pipeline.coarse.mobilePositionEstimate, cfg);
result.refinedErrors = position_errors( ...
    pipeline.refined.staticPositionEstimate, ...
    pipeline.refined.mobilePositionEstimate, cfg);
result.residualComponents = residual_components( ...
    pipeline.refined.finalResidual);
result.structuralPass = isequal(size(YSeparated), [40, 40, 3]) && ...
    isequal(size(baseReflections), [100, 40, 3]) && ...
    pipeline.coarse.numGridCandidates == 1331 && ...
    all(isfinite([pipeline.channel.delayEstimate; ...
    pipeline.channel.xiEstimate; pipeline.channel.zetaEstimate; ...
    pipeline.refined.staticPositionEstimate; ...
    pipeline.refined.mobilePositionEstimate])) && ...
    isfinite(result.reconstructionRelativeResidual);
end

function result = validate_small_noisy(cfgTemplate, resultsRoot)
snrValues = [0, 30];
seeds = (20260918:20260922).';
records = repmat(empty_record(), numel(seeds) * numel(snrValues), 1);
observations = cell(numel(seeds), numel(snrValues));
baseProfiles = cell(numel(seeds), numel(snrValues));
standardizedNoise = cell(numel(seeds), numel(snrValues));
recordIndex = 0;

for snrIndex = 1:numel(snrValues)
    for seedIndex = 1:numel(seeds)
        recordIndex = recordIndex + 1;
        seed = seeds(seedIndex);
        cfg = cfgTemplate;
        cfg.noise.snr_db = snrValues(snrIndex);
        rng(seed, 'twister');
        timer = tic;
        current = empty_record();
        current.seed = seed;
        current.snrDb = snrValues(snrIndex);
        trialFile = fullfile(resultsRoot, sprintf( ...
            'trial_seed%d_snr%+03ddB.mat', seed, snrValues(snrIndex)));
        try
            [YSeparated, noiseVariance, baseReflections] = ...
                generate_noisy_separated_signal(cfg);
            [~, truth] = generate_noiseless_signal_without_rng_change(cfg);
            noiselessSeparated = synthesize_from_parameters( ...
                truth.tau, truth.xi, truth.zeta, truth.rho, ...
                baseReflections, (1:cfg.ris.count).', cfg);
            zNoise = (YSeparated - noiselessSeparated) / sqrt(noiseVariance);
            observations{seedIndex, snrIndex} = YSeparated;
            baseProfiles{seedIndex, snrIndex} = baseReflections;
            standardizedNoise{seedIndex, snrIndex} = zNoise;

            pipeline = run_paper_baseline_A( ...
                YSeparated, noiseVariance, baseReflections, cfg);
            errors = position_errors(pipeline.refined.staticPositionEstimate, ...
                pipeline.refined.mobilePositionEstimate, cfg);
            coarseErrors = position_errors( ...
                pipeline.coarse.staticPositionEstimate, ...
                pipeline.coarse.mobilePositionEstimate, cfg);
            finiteOutput = all(isfinite([ ...
                pipeline.refined.staticPositionEstimate; ...
                pipeline.refined.mobilePositionEstimate]));
            normalExit = all(pipeline.channel.exitflag > 0) && ...
                pipeline.refined.exitflag > 0;
            if ~finiteOutput
                current.category = "1_no_valid_output";
            elseif ~normalExit
                current.category = "2_finite_optimizer_not_normal";
            else
                current.category = ...
                    "3_finite_normal_return_error_recorded";
            end
            current.hasFiniteOutput = finiteOutput;
            current.normalExit = normalExit;
            current.numValidCandidates = pipeline.coarse.numValidCandidates;
            current.ue1CoarseError = coarseErrors.ue1;
            current.ue2CoarseError = coarseErrors.ue2;
            current.jointCoarseError = coarseErrors.joint;
            current.ue1Error = errors.ue1;
            current.ue2Error = errors.ue2;
            current.jointError = errors.joint;
            current.channelExitflag = mat2str(pipeline.channel.exitflag.');
            current.channelIterations = mat2str(pipeline.channel.iterations.');
            current.channelInitialCost = ...
                sum(pipeline.channel.initialResidualCost);
            current.channelFinalCost = ...
                sum(pipeline.channel.finalResidualCost);
            current.positionExitflag = pipeline.refined.exitflag;
            current.positionIterations = pipeline.refined.iterations;
            current.positionInitialCost = pipeline.refined.initialCost;
            current.positionFinalCost = pipeline.refined.finalCost;
            current.positionFirstOrderOptimality = ...
                pipeline.refined.firstOrderOptimality;
            components = residual_components(pipeline.refined.finalResidual);
            current.finalXiCost = components.xiSquared;
            current.finalZetaCost = components.zetaSquared;
            current.finalTauCost = components.tauSquared;
            current.exceptionIdentifier = "";
            current.exceptionMessage = "";
            generationInfo.signalPower = ...
                noiseVariance * cfg.dft.num_blocks * 10^(cfg.noise.snr_db/10);
            generationInfo.separatedNoiseVariance = noiseVariance;
            generationInfo.standardizedSeparatedNoise = zNoise;
            save(trialFile, 'cfg', 'seed', 'YSeparated', 'noiseVariance', ...
                'baseReflections', 'generationInfo', 'pipeline', ...
                'errors', 'coarseErrors', '-v7.3');
        catch exception
            current.category = "1_no_valid_output";
            current.exceptionIdentifier = string(exception.identifier);
            current.exceptionMessage = string(exception.message);
            save(trialFile, 'cfg', 'seed', 'exception', '-v7.3');
        end
        current.elapsedSeconds = toc(timer);
        save(trialFile, 'current', '-append');
        records(recordIndex) = current;
    end
end

profileDifferences = zeros(numel(seeds), 1);
noiseDifferences = zeros(numel(seeds), 1);
for seedIndex = 1:numel(seeds)
    if ~isempty(baseProfiles{seedIndex, 1}) && ...
            ~isempty(baseProfiles{seedIndex, 2})
        profileDifferences(seedIndex) = max(abs( ...
            baseProfiles{seedIndex, 1}(:) - baseProfiles{seedIndex, 2}(:)));
        noiseDifferences(seedIndex) = max(abs( ...
            standardizedNoise{seedIndex, 1}(:) - ...
            standardizedNoise{seedIndex, 2}(:)));
    else
        profileDifferences(seedIndex) = nan;
        noiseDifferences(seedIndex) = nan;
    end
end

result.snrDbValues = snrValues;
result.seeds = seeds;
result.summary = records;
result.maxBaseReflectionDifferenceAcrossSnr = profileDifferences;
result.maxStandardizedNoiseDifferenceAcrossSnr = noiseDifferences;
result.sharedRandomnessVerified = all(profileDifferences == 0) && ...
    all(noiseDifferences < 1e-12);
result.observationSizes = cellfun(@size_as_text, observations, ...
    'UniformOutput', false);
end

function [Y, truth] = generate_noiseless_signal_without_rng_change(cfg)
state = rng;
cleanup = onCleanup(@() rng(state)); %#ok<NASGU>
[Y, truth] = generate_noiseless_signal(cfg);
end

function output = reconstruct_separated(channel, baseReflections, active, cfg)
output = synthesize_from_parameters(channel.delayEstimate, ...
    channel.xiEstimate, channel.zetaEstimate, channel.rhoEstimate, ...
    baseReflections, active, cfg);
end

function output = synthesize_from_parameters( ...
    tau, xi, zeta, rho, baseReflections, active, cfg)
output = complex(zeros(cfg.ofdm.num_subcarriers, ...
    size(baseReflections, 2), cfg.ris.count));
n = (0:(cfg.ofdm.num_subcarriers - 1)).';
for localIndex = 1:numel(active)
    ell = active(localIndex);
    delay = exp(-1i * 2 * pi * n * cfg.ofdm.subcarrier_spacing * ...
        tau(localIndex));
    spatial = baseReflections(:, :, ell).' * ...
        upa_steering(xi(localIndex), zeta(localIndex), cfg);
    output(:, :, ell) = rho(localIndex) * delay * spatial.';
end
end

function value = finite_difference_jacobian_check( ...
    staticPosition, mobilePosition, tau, xi, zeta, active, cfg, analytic)
x = [staticPosition; mobilePosition];
step = 1e-4;
numeric = zeros(size(analytic));
for column = 1:numel(x)
    xp = x; xm = x;
    xp(column) = xp(column) + step;
    xm(column) = xm(column) - step;
    numeric(:, column) = (parameter_residual(xp, tau, xi, zeta, ...
        active, cfg) - parameter_residual(xm, tau, xi, zeta, ...
        active, cfg)) / (2 * step);
end
angleRows = reshape([3*(0:numel(active)-1)+1; ...
    3*(0:numel(active)-1)+2], [], 1);
delayRows = 3*(1:numel(active));
angleDifference = analytic(angleRows, :) - numeric(angleRows, :);
delayDifference = analytic(delayRows, :) - numeric(delayRows, :);
value.stepMetres = step;
value.angleMaxAbsError = max(abs(angleDifference), [], 'all');
value.delayMaxAbsError = max(abs(delayDifference), [], 'all');
value.angleRelativeFroError = norm(angleDifference, 'fro') / ...
    max(norm(numeric(angleRows, :), 'fro'), eps);
value.delayRelativeFroError = norm(delayDifference, 'fro') / ...
    max(norm(numeric(delayRows, :), 'fro'), eps);
value.analytic = analytic;
value.numeric = numeric;
end

function residual = parameter_residual(x, tau, xi, zeta, active, cfg)
residual = zeros(3 * numel(active), 1);
for index = 1:numel(active)
    ell = active(index);
    ds = x(1:3) - cfg.geometry.ris_centers(:, ell);
    dm = x(4:6) - cfg.geometry.ris_centers(:, ell);
    us = ds / norm(ds);
    um = dm / norm(dm);
    predictedXi = cfg.geometry.ris_axis_x(:, ell).' * (us + um);
    predictedZeta = cfg.geometry.ris_axis_y(:, ell).' * (us + um);
    predictedTau = (norm(ds) + norm(dm)) / cfg.constants.c;
    rows = 3 * (index - 1) + (1:3);
    residual(rows) = [xi(index) - predictedXi; ...
        zeta(index) - predictedZeta; tau(index) - predictedTau];
end
end

function value = channel_errors(channel, truth, active)
value.delayDirect = channel.delayEstimate - truth.tau(active);
value.delayPeriodic = periodic_difference(channel.delayEstimate, ...
    truth.tau(active), 1 / 2.5e6);
value.xiDirect = channel.xiEstimate - truth.xi(active);
value.xiPeriodic = periodic_difference(channel.xiEstimate, ...
    truth.xi(active), 2);
value.zetaDirect = channel.zetaEstimate - truth.zeta(active);
value.zetaPeriodic = periodic_difference(channel.zetaEstimate, ...
    truth.zeta(active), 2);
end

function value = position_errors(p1, p2, cfg)
value.ue1 = norm(p1 - cfg.geometry.static_user);
value.ue2 = norm(p2 - cfg.geometry.mobile_user);
value.joint = hypot(value.ue1, value.ue2);
end

function value = residual_components(residual)
value.xi = residual(1:3:end);
value.zeta = residual(2:3:end);
value.tauSeconds = residual(3:3:end);
value.xiSquared = sum(value.xi.^2);
value.zetaSquared = sum(value.zeta.^2);
value.tauSquared = sum(value.tauSeconds.^2);
value.totalSquared = sum(residual.^2);
end

function value = normal_projections(position, active, cfg)
value = zeros(numel(active), 1);
for index = 1:numel(active)
    ell = active(index);
    direction = position - cfg.geometry.ris_centers(:, ell);
    direction = direction / norm(direction);
    value(index) = cfg.geometry.ris_axis_normal(:, ell).' * direction;
end
end

function value = periodic_wrap(value, period)
value = mod(value + period / 2, period) - period / 2;
end

function value = periodic_difference(first, second, period)
value = mod(first - second + period / 2, period) - period / 2;
end

function value = search_info(cfg)
value.lowerBound = cfg.paperBaselineA.lowerBound;
value.upperBound = cfg.paperBaselineA.upperBound;
value.gridStep = cfg.paperBaselineA.gridStep;
end

function value = empty_record()
value = struct('seed', nan, 'snrDb', nan, 'category', "", ...
    'hasFiniteOutput', false, 'normalExit', false, ...
    'numValidCandidates', nan, 'ue1CoarseError', nan, ...
    'ue2CoarseError', nan, 'jointCoarseError', nan, ...
    'ue1Error', nan, 'ue2Error', nan, 'jointError', nan, ...
    'channelExitflag', "", 'channelIterations', "", ...
    'channelInitialCost', nan, 'channelFinalCost', nan, ...
    'positionExitflag', nan, 'positionIterations', nan, ...
    'positionInitialCost', nan, 'positionFinalCost', nan, ...
    'positionFirstOrderOptimality', nan, 'finalXiCost', nan, ...
    'finalZetaCost', nan, 'finalTauCost', nan, ...
    'exceptionIdentifier', "", 'exceptionMessage', "", ...
    'elapsedSeconds', nan);
end

function value = size_as_text(x)
if isempty(x), value = ''; else, value = mat2str(size(x)); end
end

function write_text_summary(validation, fileName)
fid = fopen(fileName, 'w');
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '%s\n', validation.label);
fprintf(fid, 'Created: %s\nMATLAB: %s\n\n', ...
    validation.createdAt, validation.matlabVersion);
fprintf(fid, 'A structural pass: %d\n', validation.A.structuralPass);
fprintf(fid, 'A candidates valid/total: %d/%d\n', ...
    validation.A.coarse.numValidCandidates, ...
    validation.A.coarse.numGridCandidates);
fprintf(fid, 'A coarse errors UE1/UE2/joint: %.12g %.12g %.12g\n', ...
    validation.A.coarseErrors.ue1, validation.A.coarseErrors.ue2, ...
    validation.A.coarseErrors.joint);
fprintf(fid, 'A refined errors UE1/UE2/joint: %.12g %.12g %.12g\n', ...
    validation.A.refinedErrors.ue1, validation.A.refinedErrors.ue2, ...
    validation.A.refinedErrors.joint);
fprintf(fid, 'A Jacobian angle/delay max abs: %.12g %.12g\n\n', ...
    validation.A.jacobianCheck.angleMaxAbsError, ...
    validation.A.jacobianCheck.delayMaxAbsError);
fprintf(fid, 'B structural pass: %d\n', validation.B.structuralPass);
fprintf(fid, 'B reconstruction relative residual: %.12g\n', ...
    validation.B.reconstructionRelativeResidual);
fprintf(fid, 'B refined errors UE1/UE2/joint: %.12g %.12g %.12g\n\n', ...
    validation.B.refinedErrors.ue1, validation.B.refinedErrors.ue2, ...
    validation.B.refinedErrors.joint);
fprintf(fid, 'C shared randomness verified: %d\n', ...
    validation.C.sharedRandomnessVerified);
for index = 1:numel(validation.C.summary)
    r = validation.C.summary(index);
    fprintf(fid, ['C seed=%d SNR=%g category=%s valid=%g ', ...
        'errors=%.12g/%.12g/%.12g exit=%g elapsed=%.6g exception=%s\n'], ...
        r.seed, r.snrDb, r.category, r.numValidCandidates, ...
        r.ue1Error, r.ue2Error, r.jointError, r.positionExitflag, ...
        r.elapsedSeconds, r.exceptionIdentifier);
end
end

function restore_state(originalPath, originalRng)
path(originalPath);
rng(originalRng);
end
