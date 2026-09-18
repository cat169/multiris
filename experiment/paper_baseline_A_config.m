function cfg = paper_baseline_A_config()
%PAPER_BASELINE_A_CONFIG Historical physical model with authorized TTT support.

cfg = simulation_config();

expectedCenters = [0, 0, -15; 20, -20, 5; 10, 10, 4.5];
assert(isequal(cfg.geometry.static_user, [-8; -6; 1.5]));
assert(isequal(cfg.geometry.mobile_user, [7; 7; 1.5]));
assert(isequal(cfg.geometry.ris_centers, expectedCenters));
assert(cfg.ris.count == 3 && cfg.ris.nx == 10 && cfg.ris.ny == 10);
assert(cfg.ris.num_elements == 100);
assert(cfg.constants.lambda == 10.7e-3);
assert(cfg.ris.dx == cfg.constants.lambda / 2);
assert(cfg.ris.dy == cfg.constants.lambda / 2);
assert(cfg.ofdm.num_subcarriers == 40);
assert(cfg.ofdm.bandwidth == 100e6);
assert(cfg.ofdm.subcarrier_spacing == 100e6 / 40);
assert(cfg.dft.num_blocks == 3);
assert(cfg.dft.base_configs_per_block == 40);
assert(cfg.dft.num_symbols == 120);

% Authorized experiment-only override. The historical snapshot is TFT.
cfg.paths.visible = [true; true; true];

% Truth-assisted coarse region for comparison A. It is consumed only by the
% experiment driver and is not passed into channel estimation.
cfg.paperBaselineA.lowerBound = cfg.geometry.static_user - ones(3, 1);
cfg.paperBaselineA.upperBound = cfg.geometry.static_user + ones(3, 1);
cfg.paperBaselineA.gridStep = 0.2;
cfg.paperBaselineA.oracleSupport = (1:cfg.ris.count).';
cfg.paperBaselineA.label = [ ...
    "Multi-RIS 3D adapted to the single-frame two-user scenario, " + ...
    "with oracle support and a truth-centered coarse search region."];
end
