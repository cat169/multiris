function result = refine_position_map( ...
    channelEstimates, channelFimResults, positionPriorMessage, ...
    initialPosition, activeRisIndices, cfg)
%REFINE_POSITION_MAP Fisher-weighted six-dimensional position MAP estimate.
% Residual order for every RIS is [xi; zeta; nu], nu = Delta_f * tau.

activeRisIndices = activeRisIndices(:);
priorMean = positionPriorMessage.mean(:);
priorCovariance = positionPriorMessage.covariance;
if isempty(initialPosition)
    x0 = priorMean;
else
    x0 = initialPosition(:);
end

numActiveRis = numel(activeRisIndices);
assert(numActiveRis > 0 && numel(priorMean) == 6 && numel(x0) == 6);
assert(isequal(size(priorCovariance), [6, 6]));
assert(numel(channelEstimates.xiEstimate) == numActiveRis && ...
    numel(channelEstimates.zetaEstimate) == numActiveRis && ...
    numel(channelEstimates.delayEstimate) == numActiveRis);
assert(numel(channelFimResults) == numActiveRis);

% C0 = L0*L0.'; ||L0\(x-m0)||^2 = (x-m0).'*C0^(-1)*(x-m0).
L0 = chol((priorCovariance + priorCovariance.') / 2, 'lower');
priorSquareRoot = L0 \ eye(6);
fisherSquareRoots = cell(numActiveRis, 1);
for activeIndex = 1:numActiveRis
    if iscell(channelFimResults)
        fisherResult = channelFimResults{activeIndex};
    else
        fisherResult = channelFimResults(activeIndex);
    end
    fisher = (fisherResult.fim + fisherResult.fim.') / 2;
    [eigenvectors, eigenvalueMatrix] = eig(fisher);
    eigenvalues = real(diag(eigenvalueMatrix));
    tolerance = 100 * eps * max(norm(fisher, 2), realmin);
    if any(eigenvalues < -tolerance)
        error('refine_position_map:IndefiniteFisher', ...
            'A channel Fisher information matrix is not positive semidefinite.');
    end
    eigenvalues(eigenvalues < 0) = 0;
    % R.'*R = F, including zero-information eigendirections.
    fisherSquareRoots{activeIndex} = ...
        diag(sqrt(eigenvalues)) * eigenvectors.';
end

options = optimoptions('lsqnonlin', ...
    'Algorithm', 'trust-region-reflective', ...
    'SpecifyObjectiveGradient', true, ...
    'Display', 'off', ...
    'MaxIterations', 500, ...
    'MaxFunctionEvaluations', 5000, ...
    'FunctionTolerance', 1e-12, ...
    'StepTolerance', 1e-12, ...
    'OptimalityTolerance', 1e-10);

[positionEstimate, resnorm, ~, exitflag, output] = lsqnonlin( ...
    @mapResidual, x0, [], [], options);

result.positionEstimate = positionEstimate;
result.staticPositionEstimate = positionEstimate(1:3);
result.mobilePositionEstimate = positionEstimate(4:6);
result.objective = 0.5 * resnorm;
result.exitflag = exitflag;
result.iterations = output.iterations;

    function [residual, jacobian] = mapResidual(x)
        [residual, jacobian] = positionMapResidual( ...
            x, channelEstimates, fisherSquareRoots, priorMean, ...
            priorSquareRoot, activeRisIndices, cfg);
    end
end

function [residual, jacobian] = positionMapResidual( ...
    x, channelEstimates, fisherSquareRoots, priorMean, ...
    priorSquareRoot, activeRisIndices, cfg)
numActiveRis = numel(activeRisIndices);
residual = zeros(3 * numActiveRis + 6, 1);
jacobian = zeros(3 * numActiveRis + 6, 6);
identity3 = eye(3);
staticPosition = x(1:3);
mobilePosition = x(4:6);
periods = [2; 2; 1];
branchOffsets = allBranchOffsets();

for activeIndex = 1:numActiveRis
    risIndex = activeRisIndices(activeIndex);
    risCenter = cfg.geometry.ris_centers(:, risIndex);
    staticDisplacement = staticPosition - risCenter;
    mobileDisplacement = mobilePosition - risCenter;
    staticDistance = norm(staticDisplacement);
    mobileDistance = norm(mobileDisplacement);
    if staticDistance == 0 || mobileDistance == 0
        error('refine_position_map:PositionAtRisCenter', ...
            'A candidate position coincides with an active RIS center.');
    end

    staticDirection = staticDisplacement / staticDistance;
    mobileDirection = mobileDisplacement / mobileDistance;
    axisXi = cfg.geometry.ris_axis_x(:, risIndex);
    axisZeta = cfg.geometry.ris_axis_y(:, risIndex);
    directionSum = staticDirection + mobileDirection;
    predictedParameters = [axisXi.' * directionSum; ...
        axisZeta.' * directionSum; ...
        cfg.ofdm.subcarrier_spacing * ...
        (staticDistance + mobileDistance) / cfg.constants.c];
    estimatedParameters = [channelEstimates.xiEstimate(activeIndex); ...
        channelEstimates.zetaEstimate(activeIndex); ...
        cfg.ofdm.subcarrier_spacing * ...
        channelEstimates.delayEstimate(activeIndex)];

    rawResidual = estimatedParameters - predictedParameters;
    centerBranch = -round(rawResidual ./ periods);
    candidateResiduals = rawResidual + periods .* ...
        (centerBranch + branchOffsets);
    fisherSquareRoot = fisherSquareRoots{activeIndex};
    weightedCandidates = fisherSquareRoot * candidateResiduals;
    [~, bestBranch] = min(sum(weightedCandidates.^2, 1));

    staticDirectionJacobian = ...
        (identity3 - staticDirection * staticDirection.') / staticDistance;
    mobileDirectionJacobian = ...
        (identity3 - mobileDirection * mobileDirection.') / mobileDistance;
    predictedJacobian = [ ...
        axisXi.' * staticDirectionJacobian, ...
        axisXi.' * mobileDirectionJacobian; ...
        axisZeta.' * staticDirectionJacobian, ...
        axisZeta.' * mobileDirectionJacobian; ...
        cfg.ofdm.subcarrier_spacing * staticDirection.' / cfg.constants.c, ...
        cfg.ofdm.subcarrier_spacing * mobileDirection.' / cfg.constants.c];

    rows = 3 * (activeIndex - 1) + (1:3);
    residual(rows) = weightedCandidates(:, bestBranch);
    jacobian(rows, :) = -fisherSquareRoot * predictedJacobian;
end

priorRows = 3 * numActiveRis + (1:6);
residual(priorRows) = priorSquareRoot * (x - priorMean);
jacobian(priorRows, :) = priorSquareRoot;
end

function offsets = allBranchOffsets()
% Columns enumerate {-1,0,1}^3 around the coordinate-wise nearest branch.
[offsetXi, offsetZeta, offsetNu] = ndgrid(-1:1, -1:1, -1:1);
offsets = [offsetXi(:), offsetZeta(:), offsetNu(:)].';
end
