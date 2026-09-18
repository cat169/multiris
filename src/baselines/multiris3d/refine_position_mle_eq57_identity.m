function result = refine_position_mle_eq57_identity( ...
    delayEstimate, xiEstimate, zetaEstimate, activeRisIndices, ...
    initialStaticPosition, initialMobilePosition, cfg)
%REFINE_POSITION_MLE_EQ57_IDENTITY Eq. (57) with Sigma_eta equal to I.
% The per-RIS residual is [xiHat-xiPred; zetaHat-zetaPred;
% tauHat-tauPred], where tau is expressed in seconds. The overall sign is
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
[initialResidual, initialJacobian] = positionResidualEq57Identity( ...
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
result.residualOrder = ["xi"; "zeta"; "tau_seconds"];

    function [residual, jacobian] = residualForSolver(x)
        [residual, jacobian] = positionResidualEq57Identity( ...
            x, delayEstimate, xiEstimate, zetaEstimate, ...
            activeRisIndices, cfg);
    end
end

function [residual, jacobian] = positionResidualEq57Identity( ...
    x, delayEstimate, xiEstimate, zetaEstimate, activeRisIndices, cfg)
numActiveRis = numel(activeRisIndices);
residual = zeros(3 * numActiveRis, 1);
jacobian = zeros(3 * numActiveRis, 6);
identity3 = eye(3);
staticPosition = x(1:3);
mobilePosition = x(4:6);

for activeIndex = 1:numActiveRis
    risIndex = activeRisIndices(activeIndex);
    risCenter = cfg.geometry.ris_centers(:, risIndex);
    staticDisplacement = staticPosition - risCenter;
    mobileDisplacement = mobilePosition - risCenter;
    staticDistance = norm(staticDisplacement);
    mobileDistance = norm(mobileDisplacement);
    if staticDistance == 0 || mobileDistance == 0
        error('refine_position_mle_eq57_identity:PositionAtRisCenter', ...
            'A candidate position coincides with an active RIS center.');
    end

    staticDirection = staticDisplacement / staticDistance;
    mobileDirection = mobileDisplacement / mobileDistance;
    axisXi = cfg.geometry.ris_axis_x(:, risIndex);
    axisZeta = cfg.geometry.ris_axis_y(:, risIndex);
    directionSum = staticDirection + mobileDirection;
    predictedXi = axisXi.' * directionSum;
    predictedZeta = axisZeta.' * directionSum;
    predictedDelay = ...
        (staticDistance + mobileDistance) / cfg.constants.c;
    rows = 3 * (activeIndex - 1) + (1:3);

    residual(rows) = [xiEstimate(activeIndex) - predictedXi; ...
        zetaEstimate(activeIndex) - predictedZeta; ...
        delayEstimate(activeIndex) - predictedDelay];

    staticDirectionJacobian = ...
        (identity3 - staticDirection * staticDirection.') / ...
        staticDistance;
    mobileDirectionJacobian = ...
        (identity3 - mobileDirection * mobileDirection.') / ...
        mobileDistance;
    jacobian(rows, :) = -[ ...
        axisXi.' * staticDirectionJacobian, ...
        axisXi.' * mobileDirectionJacobian; ...
        axisZeta.' * staticDirectionJacobian, ...
        axisZeta.' * mobileDirectionJacobian; ...
        staticDirection.' / cfg.constants.c, ...
        mobileDirection.' / cfg.constants.c];
end
end
