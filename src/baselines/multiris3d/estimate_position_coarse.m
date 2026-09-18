function result = estimate_position_coarse( ...
    delayEstimate, xiEstimate, zetaEstimate, ...
    activeRisIndices, positionSearchInfo, cfg)
%ESTIMATE_POSITION_COARSE 执行 Multi-RIS 3D 粗位置估计。
% 在 static UE 网格上搜索；由 xi/zeta 与固定正法向恢复 mobile direction，
% 使用 bearing-line 最近点等权融合 mobile UE，再以 B=0 delay consistency 选点。

activeRisIndices = activeRisIndices(:);
delayEstimate = delayEstimate(:);
xiEstimate = xiEstimate(:);
zetaEstimate = zetaEstimate(:);
numActiveRis = numel(activeRisIndices);
if numActiveRis < 2
    error('estimate_position_coarse:InsufficientActiveRis', ...
        'At least two active RIS paths are required.');
end
assert(all(ismember(activeRisIndices, 1:cfg.ris.count)));
assert(numel(delayEstimate) == numActiveRis);
assert(numel(xiEstimate) == numActiveRis);
assert(numel(zetaEstimate) == numActiveRis);
assert(all(isfinite(delayEstimate)));
assert(all(isfinite(xiEstimate)));
assert(all(isfinite(zetaEstimate)));
assert(isfield(cfg.geometry, 'ris_axis_normal'));
assert(isequal(size(cfg.geometry.ris_axis_normal), ...
    [3, cfg.ris.count]));

lowerBound = positionSearchInfo.lowerBound(:);
upperBound = positionSearchInfo.upperBound(:);
gridStep = positionSearchInfo.gridStep;
assert(isequal(size(lowerBound), [3, 1]));
assert(isequal(size(upperBound), [3, 1]));
assert(all(isfinite(lowerBound)));
assert(all(isfinite(upperBound)));
assert(all(lowerBound <= upperBound));
assert(isscalar(gridStep) && isfinite(gridStep) && gridStep > 0);

xGrid = (lowerBound(1):gridStep:upperBound(1)).';
yGrid = (lowerBound(2):gridStep:upperBound(2)).';
zGrid = (lowerBound(3):gridStep:upperBound(3)).';
[gridX, gridY, gridZ] = ndgrid(xGrid, yGrid, zGrid);
staticGrid = [gridX(:).'; gridY(:).'; gridZ(:).'];
numGridCandidates = size(staticGrid, 2);
candidateCost = inf(numGridCandidates, 1);
numericTolerance = 1e-12;

% 每个 static candidate 必须对全部 active RIS 和全部 RIS pair 同时可行。
for gridIndex = 1:numGridCandidates
    staticCandidate = staticGrid(:, gridIndex);
    geometry = recoverMobileGeometry( ...
        staticCandidate, xiEstimate, zetaEstimate, ...
        activeRisIndices, numericTolerance, cfg);
    if geometry.isValid
        delayResidual = calculateDelayResidual( ...
            staticCandidate, geometry.mobilePosition, ...
            delayEstimate, activeRisIndices, cfg);
        candidateCost(gridIndex) = sum(abs(delayResidual));
    end
end

numValidCandidates = sum(isfinite(candidateCost));
[minimumCost, bestGridIndex] = min(candidateCost);
if ~isfinite(minimumCost)
    error('estimate_position_coarse:NoValidCandidate', ...
        'All static grid candidates have invalid geometry.');
end

staticPositionEstimate = staticGrid(:, bestGridIndex);
bestGeometry = recoverMobileGeometry( ...
    staticPositionEstimate, xiEstimate, zetaEstimate, ...
    activeRisIndices, numericTolerance, cfg);
assert(bestGeometry.isValid);
bestDelayResidual = calculateDelayResidual( ...
    staticPositionEstimate, bestGeometry.mobilePosition, ...
    delayEstimate, activeRisIndices, cfg);

result.staticPositionEstimate = staticPositionEstimate;
result.mobilePositionEstimate = bestGeometry.mobilePosition;
result.minimumCost = minimumCost;
result.bestGridIndex = bestGridIndex;
result.bestDelayResidual = bestDelayResidual;
result.searchLowerBound = lowerBound;
result.searchUpperBound = upperBound;
result.gridStep = gridStep;
result.gridSize = [numel(xGrid); numel(yGrid); numel(zGrid)];
result.numGridCandidates = numGridCandidates;
result.numValidCandidates = numValidCandidates;
result.activeRisIndices = activeRisIndices;
result.bestRecoveredMobileDirections = bestGeometry.mobileDirections;
result.bestPairRisIndices = bestGeometry.pairRisIndices;
result.bestPairPositionEstimates = bestGeometry.pairPositionEstimates;
result.bestPairRanges = bestGeometry.pairRanges;
result.bestPairClosestDistances = bestGeometry.pairClosestDistances;
result.bestPairConditionNumbers = bestGeometry.pairConditionNumbers;
end

function geometry = recoverMobileGeometry( ...
    staticCandidate, xiEstimate, zetaEstimate, ...
    activeRisIndices, numericTolerance, cfg)
numActiveRis = numel(activeRisIndices);
mobileDirections = zeros(3, numActiveRis);
geometry = emptyGeometry(numActiveRis, activeRisIndices);

% xi/zeta 给出面内分量，第三分量始终沿部署配置的固定正法向取正根。
for activeIndex = 1:numActiveRis
    ell = activeRisIndices(activeIndex);
    risCenter = cfg.geometry.ris_centers(:, ell);
    deltaStatic = staticCandidate - risCenter;
    staticRange = norm(deltaStatic);
    if staticRange == 0
        return
    end
    uStatic = deltaStatic / staticRange;
    eXi = cfg.geometry.ris_axis_x(:, ell);
    eZeta = cfg.geometry.ris_axis_y(:, ell);
    eNormal = cfg.geometry.ris_axis_normal(:, ell);
    a = xiEstimate(activeIndex) - eXi' * uStatic;
    b = zetaEstimate(activeIndex) - eZeta' * uStatic;
    rPerpSquared = 1 - a^2 - b^2;
    if rPerpSquared < -numericTolerance
        return
    elseif rPerpSquared < 0
        rPerpSquared = 0;
    end
    cNormal = sqrt(rPerpSquared);
    mobileDirections(:, activeIndex) = ...
        a * eXi + b * eZeta + cNormal * eNormal;
end

% 对每个 unordered pair 使用 line i 上的最近点 pIJ，不使用 midpoint。
pairLocalIndices = nchoosek(1:numActiveRis, 2);
pairRisIndices = activeRisIndices(pairLocalIndices);
numPairs = size(pairLocalIndices, 1);
pairPositionEstimates = zeros(3, numPairs);
pairRanges = zeros(numPairs, 2);
pairClosestDistances = zeros(numPairs, 1);
pairConditionNumbers = zeros(numPairs, 1);
for pairIndex = 1:numPairs
    localI = pairLocalIndices(pairIndex, 1);
    localJ = pairLocalIndices(pairIndex, 2);
    ellI = activeRisIndices(localI);
    ellJ = activeRisIndices(localJ);
    pI = cfg.geometry.ris_centers(:, ellI);
    pJ = cfg.geometry.ris_centers(:, ellJ);
    uI = mobileDirections(:, localI);
    uJ = mobileDirections(:, localJ);
    A = [uI, -uJ];
    if rank(A) < 2
        return
    end
    ranges = A \ (pJ - pI);
    pointI = pI + ranges(1) * uI;
    pointJ = pJ + ranges(2) * uJ;
    pairPositionEstimates(:, pairIndex) = pointI;
    pairRanges(pairIndex, :) = ranges.';
    pairClosestDistances(pairIndex) = norm(pointI - pointJ);
    pairConditionNumbers(pairIndex) = cond(A' * A);
end

% Eq. (53) 在无额外权重时对所有 pair 采用等权 receiver fusion。
geometry.isValid = true;
geometry.mobileDirections = mobileDirections;
geometry.pairRisIndices = pairRisIndices;
geometry.pairPositionEstimates = pairPositionEstimates;
geometry.pairRanges = pairRanges;
geometry.pairClosestDistances = pairClosestDistances;
geometry.pairConditionNumbers = pairConditionNumbers;
geometry.mobilePosition = mean(pairPositionEstimates, 2);
end

function geometry = emptyGeometry(numActiveRis, activeRisIndices)
numPairs = numActiveRis * (numActiveRis - 1) / 2;
geometry.isValid = false;
geometry.mobileDirections = nan(3, numActiveRis);
geometry.pairRisIndices = nan(numPairs, 2);
geometry.pairPositionEstimates = nan(3, numPairs);
geometry.pairRanges = nan(numPairs, 2);
geometry.pairClosestDistances = nan(numPairs, 1);
geometry.pairConditionNumbers = nan(numPairs, 1);
geometry.mobilePosition = nan(3, 1);
if numPairs > 0
    pairLocalIndices = nchoosek(1:numActiveRis, 2);
    geometry.pairRisIndices = activeRisIndices(pairLocalIndices);
end
end

function delayResidual = calculateDelayResidual( ...
    staticPosition, mobilePosition, ...
    delayEstimate, activeRisIndices, cfg)
risCenters = cfg.geometry.ris_centers(:, activeRisIndices);
staticDistance = vecnorm(staticPosition - risCenters, 2, 1).';
mobileDistance = vecnorm(mobilePosition - risCenters, 2, 1).';
% 当前 shared 场景无 LOS 且 B=0，所有 active RIS 使用相同 delay 权重。
delayResidual = staticDistance + mobileDistance - ...
    cfg.constants.c * delayEstimate;
end
