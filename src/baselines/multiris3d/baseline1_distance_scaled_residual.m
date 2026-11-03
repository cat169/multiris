function [residual, jacobian] = baseline1_distance_scaled_residual( ...
    x, delayEstimate, xiEstimate, zetaEstimate, activeRisIndices, cfg)
%BASELINE1_DISTANCE_SCALED_RESIDUAL Public pure evaluator for the C variant.
% Geometry and residual/Jacobian only, no optimizer or truth-based selection.
x=x(:);delayEstimate=delayEstimate(:);xiEstimate=xiEstimate(:);zetaEstimate=zetaEstimate(:);
activeRisIndices=activeRisIndices(:);
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
        error('refine_position_mle_eq57_distance_scaled:PositionAtRisCenter', ...
            'A candidate position coincides with an active RIS center.');
    end

    staticDirection = staticDisplacement / staticDistance;
    mobileDirection = mobileDisplacement / mobileDistance;
    axisXi = cfg.geometry.ris_axis_x(:, risIndex);
    axisZeta = cfg.geometry.ris_axis_y(:, risIndex);
    directionSum = staticDirection + mobileDirection;
    predictedXi = axisXi.' * directionSum;
    predictedZeta = axisZeta.' * directionSum;
    rows = 3 * (activeIndex - 1) + (1:3);

    residual(rows) = [xiEstimate(activeIndex) - predictedXi; ...
        zetaEstimate(activeIndex) - predictedZeta; ...
        cfg.constants.c * delayEstimate(activeIndex) - staticDistance - mobileDistance];

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
        staticDirection.', ...
        mobileDirection.'];
end
end
