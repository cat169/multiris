function [residual, jacobian] = b2_s4_residual( ...
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
