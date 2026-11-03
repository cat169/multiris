function result = refine_position_mle_eq57_distance_scaled( ...
    delayEstimate, xiEstimate, zetaEstimate, activeRisIndices, ...
    initialStaticPosition, initialMobilePosition, cfg)
%REFINE_POSITION_MLE_EQ57_DISTANCE_SCALED Independent distance-scaled Eq. (57) experiment.
% The per-RIS residual is [xiHat-xiPred; zetaHat-zetaPred;
% c*tauHat-dT-dR], where the delay residual is in metres. The overall sign is
% immaterial to least squares and matches the formal Step 8 convention.

delayEstimate = delayEstimate(:);
xiEstimate = xiEstimate(:);
zetaEstimate = zetaEstimate(:);
activeRisIndices = activeRisIndices(:);
initialStaticPosition = initialStaticPosition(:);
initialMobilePosition = initialMobilePosition(:);

numActiveRis = numel(activeRisIndices);
assert(numel(delayEstimate) == numActiveRis);
assert(numel(xiEstimate) == numActiveRis);
assert(numel(zetaEstimate) == numActiveRis);
assert(numActiveRis > 0);
assert(all(ismember(activeRisIndices, 1:cfg.ris.count)));
assert(numel(initialStaticPosition) == 3);
assert(numel(initialMobilePosition) == 3);
assert(all(isfinite(delayEstimate)));
assert(all(isfinite(xiEstimate)));
assert(all(isfinite(zetaEstimate)));
assert(all(isfinite(initialStaticPosition)));
assert(all(isfinite(initialMobilePosition)));

x0 = [initialStaticPosition; initialMobilePosition];
[initialResidual, initialJacobian] = baseline1_distance_scaled_residual( ...
    x0, delayEstimate, xiEstimate, zetaEstimate, ...
    activeRisIndices, cfg);

options = optimoptions('lsqnonlin', ...
    'Algorithm', 'trust-region-reflective', ...
    'SpecifyObjectiveGradient', true, ...
    'Display', 'off', ...
    'MaxIterations', 500, ...
    'MaxFunctionEvaluations', 5000, ...
    'FunctionTolerance', 1e-12, ...
    'StepTolerance', 1e-12, ...
    'OptimalityTolerance', 1e-10);

[xOpt, resnorm, finalResidual, exitflag, output] = lsqnonlin( ...
    @residualForSolver, x0, [], [], options);

result.staticPositionEstimate = xOpt(1:3);
result.mobilePositionEstimate = xOpt(4:6);
result.initialCost = sum(initialResidual.^2);
result.finalCost = resnorm;
result.exitflag = exitflag;
result.iterations = output.iterations;
result.funcCount = output.funcCount;
result.firstOrderOptimality = output.firstorderopt;
result.initialJacobian = initialJacobian;
result.initialResidual = initialResidual;
result.finalResidual = finalResidual;
result.residualOrder = ["xi"; "zeta"; "distance_metres"];

    function [residual, jacobian] = residualForSolver(x)
        [residual, jacobian] = baseline1_distance_scaled_residual( ...
            x, delayEstimate, xiEstimate, zetaEstimate, ...
            activeRisIndices, cfg);
    end
end
