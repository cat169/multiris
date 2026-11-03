# Baseline 2-S4 equivalence and return-contract verification

Production paths (all within one MultiRIS checkout):

```matlab
addpath(fullfile(repo,'baseline2','inference'), ...
    fullfile(repo,'src','baselines','multiris3d'), ...
    fullfile(repo,'src','structured','current','signal'));
r = infer_baseline2_frame_s4(Y,nv,W,prior,activeRisIndices,cfg);
```

`prior.mean` is 6-by-1, `prior.covariance` is 6-by-6 positive definite.
No Tracking path is required by production inference. Optimization Toolbox
is required by the unchanged original `refine_position_map`.

The retained-data check requires the original external Tracking archive and
its SHA-verified analysis functions; those are test references, not production
dependencies. The script audits 1600 saved records without optimizing them,
then recomputes 13 fixed representatives on their original seeds and compares
the old 14-candidate evaluator with the new interface. No VB/B1/PEB is called.

```matlab
addpath(fullfile(repo,'tests','baseline2_s4'));
run_baseline2_s4_equivalence(trackingRoot,outputDirectory);
results = runtests(fullfile(repo,'tests','baseline2_s4','Baseline2S4ContractTest.m'));
assert(all([results.Passed]));
```

Fixed comparison tolerances: initial positions and Fisher matrices exact;
positions 1e-8 m in infinity norm; objectives 1e-10*max(1,abs(reference));
gradient differences 1e-8*max(1,abs(reference)); indices, flags and iteration
counts exact. The original candidate-consistency tolerance remains 1e-10.

The eight class tests use isolated temporary solver doubles only for return
contracts: all valid, partial failure, all failure, original-SR failure,
failed dependent swap, finite nonconvergence, program error, invalid prior.
They do not claim numerical MAP equivalence. The numerical equivalence script
uses the real original MLE/Fisher/MAP functions and MATLAB profiler counts.

`check_baseline2_s4_standalone(repo,inputFile)` validates one frozen input in
an independent checkout with only the three production paths. `inputFile`
contains an `input` structure with Y/nv/W/prior/active/cfg, referencePosition,
referenceIndex and expectedMapCalls. Large historical data and frozen inputs
are intentionally not committed. No test runner is invoked by production.
