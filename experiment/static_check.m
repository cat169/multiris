function report = static_check()
%STATIC_CHECK Verify isolation, configuration, grids, and required toolbox.

[root, cleanup] = enter_isolated_experiment(); %#ok<ASGLU>
cfg = paper_baseline_A_config();
expected = required_function_map(root);
names = fieldnames(expected);
resolved = struct();
for index = 1:numel(names)
    name = names{index};
    resolved.(name) = which(name);
    assert(strcmpi(resolved.(name), expected.(name)), ...
        'paperA:PathLeak', '%s resolved to %s, expected %s.', ...
        name, resolved.(name), expected.(name));
end

gridCounts = round((cfg.paperBaselineA.upperBound - ...
    cfg.paperBaselineA.lowerBound) / cfg.paperBaselineA.gridStep) + 1;
assert(isequal(gridCounts, [11; 11; 11]));
assert(prod(gridCounts) == 1331);
assert(isequal(cfg.paths.visible, [true; true; true]));
assert(license('test', 'Optimization_Toolbox'), ...
    'paperA:MissingOptimizationToolbox', ...
    'Optimization Toolbox is required for lsqnonlin.');

report.root = root;
report.label = cfg.paperBaselineA.label;
report.resolvedFunctions = resolved;
report.gridCounts = gridCounts;
report.totalCandidates = prod(gridCounts);
report.visibility = cfg.paths.visible;
report.matlabVersion = version;
report.optimizationToolboxAvailable = true;
end

function [root, cleanup] = enter_isolated_experiment()
originalPath = path;
originalRng = rng;
cleanup = onCleanup(@() restore_state(originalPath, originalRng));
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
restoredefaultpath;
addpath(fullfile(root, 'src', 'structured', 'current'), '-begin');
addpath(fullfile(root, 'src', 'structured', 'current', 'signal'), '-begin');
addpath(fullfile(root, 'src', 'baselines', 'multiris3d'), '-begin');
addpath(here, '-begin');
end

function restore_state(originalPath, originalRng)
path(originalPath);
rng(originalRng);
end

function expected = required_function_map(root)
expected.simulation_config = fullfile(root, 'src', 'structured', ...
    'current', 'simulation_config.m');
expected.generate_noiseless_signal = fullfile(root, 'src', 'structured', ...
    'current', 'signal', 'generate_noiseless_signal.m');
expected.generate_noisy_separated_signal = fullfile(root, 'src', ...
    'structured', 'current', 'signal', ...
    'generate_noisy_separated_signal.m');
expected.upa_steering = fullfile(root, 'src', 'structured', 'current', ...
    'signal', 'upa_steering.m');
expected.estimate_delay_coarse = fullfile(root, 'src', 'baselines', ...
    'multiris3d', 'estimate_delay_coarse.m');
expected.estimate_spatial_coarse = fullfile(root, 'src', 'baselines', ...
    'multiris3d', 'estimate_spatial_coarse.m');
expected.refine_channel_mle = fullfile(root, 'src', 'baselines', ...
    'multiris3d', 'refine_channel_mle.m');
expected.estimate_position_coarse = fullfile(root, 'src', 'baselines', ...
    'multiris3d', 'estimate_position_coarse.m');
expected.refine_position_mle_eq57_identity = fullfile(root, 'src', ...
    'baselines', 'multiris3d', 'refine_position_mle_eq57_identity.m');
expected.paper_baseline_A_config = fullfile(root, 'experiment', ...
    'paper_baseline_A_config.m');
expected.run_paper_baseline_A = fullfile(root, 'experiment', ...
    'run_paper_baseline_A.m');
end
