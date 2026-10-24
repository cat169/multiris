%VALIDATE_CHANNEL_FIM One lightweight, read-only check on an existing trial.
% No signal generation, estimators, position optimization, RNG or file writes.
% Run this script from any working directory. Paths are restored afterwards.
baseline2Folder = fileparts(fileparts(mfilename('fullpath')));
projectFolder = fileparts(baseline2Folder);
originalPath = path;
pathCleanup = onCleanup(@() path(originalPath));
addpath(fullfile(baseline2Folder, 'inference'));
addpath(fullfile(projectFolder, 'src', 'structured', 'current', 'signal'));
sampleFile = fullfile(projectFolder, 'results', ...
    'trial_seed20260918_snr+30dB.mat');
sample = load(sampleFile, 'cfg', 'baseReflections', 'pipeline');
cfg = sample.cfg;
channel = sample.pipeline.channel;
fprintf('Existing sample: %s\n', sampleFile);
fprintf('Separated complex noise variance: %.12g\n', ...
    channel.separatedNoiseVariance);
fprintf('Parameter order: [xi; zeta; nu], nu = df*tau\n');
fprintf('RIS | derivative relative errors (xi,zeta,nu) | eigenvalues | cond\n');
derivativeErrors = zeros(numel(channel.activeRisIndices), 3);
fimResults = cell(numel(channel.activeRisIndices), 1);
h = 1e-7;
for activeIndex = 1:numel(channel.activeRisIndices)
    ell = channel.activeRisIndices(activeIndex);
    W = sample.baseReflections(:, :, ell);
    q = [channel.xiEstimate(activeIndex); ...
        channel.zetaEstimate(activeIndex); ...
        cfg.ofdm.subcarrier_spacing * channel.delayEstimate(activeIndex)];
    [info, model] = compute_channel_fim(q, ...
        channel.rhoEstimate(activeIndex), W, ...
        channel.separatedNoiseVariance, cfg);
    for parameterIndex = 1:3
        step = zeros(3, 1);
        step(parameterIndex) = h;
        difference = (referenceAtom(q + step, W, cfg) - ...
            referenceAtom(q - step, W, cfg)) / (2 * h);
        derivativeErrors(activeIndex, parameterIndex) = ...
            norm(difference - model.derivatives(:, parameterIndex)) / ...
            norm(model.derivatives(:, parameterIndex));
    end
    assert(all(derivativeErrors(activeIndex, :) < 1e-6), ...
        'Analytic derivative mismatch.');
    assert(isequal(size(info.fim), [3, 3]) && isreal(info.fim));
    assert(norm(info.fim - info.fim.', 'fro') <= info.psdTolerance);
    assert(info.isPositiveSemidefinite, 'Significant negative eigenvalue.');
    assert(isequal(info.parameterOrder, {'xi'; 'zeta'; 'nu'}));
    fprintf('%d | %.3e %.3e %.3e | %.6e %.6e %.6e | %.6g\n', ...
        ell, derivativeErrors(activeIndex, :), info.eigenvalues, ...
        info.conditionNumber);
    fprintf('    PSD tolerance %.3e; near singular: %d\n', ...
        info.psdTolerance, info.isNearlySingular);
    fimResults{activeIndex} = info;
end
fprintf('PASS: analytic derivatives, real symmetric 3x3 FIM, PSD and order.\n');
fprintf('Condition numbers reported without regularization or inversion.\n');
clear pathCleanup

function mu = referenceAtom(q, W, cfg)
% Independent reference uses the existing Baseline 1 steering function and
% its exact M(:) construction, rather than repeating analytic derivatives.
a = upa_steering(q(1), q(2), cfg);
s = W.' * a;
k = (0:(cfg.ofdm.num_subcarriers - 1)).';
d = exp(-1i * 2 * pi * k * q(3));
M = d * s.';
mu = M(:);
end
