function result = refine_channel_mle( ...
    YSeparated, baseReflections, ...
    activeRisIndices, ...
    delayEstimate, xiEstimate, zetaEstimate, cfg)
%REFINE_CHANNEL_MLE 使用 single-start fminsearch 连续精化每条 RIS 信道参数。
% optimizer 使用 q = [xi; zeta; nu]，其中 nu = subcarrierSpacing * tau，
% 并最大化归一化投影得分。所有初值均由外部 coarse estimators 显式传入。

assert(size(YSeparated, 1) == cfg.ofdm.num_subcarriers);
assert(size(YSeparated, 3) == cfg.ris.count);
assert(size(baseReflections, 1) == cfg.ris.num_elements);
assert(size(baseReflections, 3) == cfg.ris.count);
assert(size(YSeparated, 2) == size(baseReflections, 2));
assert(~isempty(activeRisIndices));
assert(all(ismember(activeRisIndices, 1:cfg.ris.count)));
assert(numel(delayEstimate) == numel(activeRisIndices));
assert(numel(xiEstimate) == numel(activeRisIndices));
assert(numel(zetaEstimate) == numel(activeRisIndices));

activeRisIndices = activeRisIndices(:);
delayEstimate = delayEstimate(:);
xiEstimate = xiEstimate(:);
zetaEstimate = zetaEstimate(:);
numActiveRis = numel(activeRisIndices);
subcarrierSpacing = cfg.ofdm.subcarrier_spacing;
subcarrierIndex = (0:(cfg.ofdm.num_subcarriers - 1)).';

options = optimset( ...
    'Display', 'off', ...
    'MaxIter', 1000, ...
    'MaxFunEvals', 3000, ...
    'TolX', 1e-10, ...
    'TolFun', 1e-12);

delayFinal = zeros(numActiveRis, 1);
xiFinal = zeros(numActiveRis, 1);
zetaFinal = zeros(numActiveRis, 1);
rhoFinal = complex(zeros(numActiveRis, 1));
initialRho = complex(zeros(numActiveRis, 1));
initialProjectionScore = zeros(numActiveRis, 1);
finalProjectionScore = zeros(numActiveRis, 1);
initialResidualCost = zeros(numActiveRis, 1);
finalResidualCost = zeros(numActiveRis, 1);
exitflag = zeros(numActiveRis, 1);
iterations = zeros(numActiveRis, 1);
funcCount = zeros(numActiveRis, 1);
optimizerObjective = zeros(numActiveRis, 1);
optimizerParameter = zeros(numActiveRis, 3);

for activeIndex = 1:numActiveRis
    ell = activeRisIndices(activeIndex);
    Hell = YSeparated(:, :, ell);
    W = baseReflections(:, :, ell);
    y = Hell(:);
    yEnergy = real(y' * y);

    xi0 = xiEstimate(activeIndex);
    zeta0 = zetaEstimate(activeIndex);
    tau0 = delayEstimate(activeIndex);
    nu0 = subcarrierSpacing * tau0;
    q0 = [xi0; zeta0; nu0];

    initial = evaluatePhysicalParameters( ...
        y, W, xi0, zeta0, tau0, subcarrierIndex, cfg);
    objectiveFunction = @(q) negativeProjectionScore( ...
        q, y, yEnergy, W, subcarrierIndex, cfg);
    [qHat, objectiveFinal, exitflagCurrent, output] = ...
        fminsearch(objectiveFunction, q0, options);

    xiFinal(activeIndex) = mod(qHat(1) + 1, 2) - 1;
    zetaFinal(activeIndex) = mod(qHat(2) + 1, 2) - 1;
    nuFinal = mod(qHat(3), 1);
    delayFinal(activeIndex) = nuFinal / subcarrierSpacing;

    final = evaluatePhysicalParameters( ...
        y, W, xiFinal(activeIndex), zetaFinal(activeIndex), ...
        delayFinal(activeIndex), subcarrierIndex, cfg);

    rhoFinal(activeIndex) = final.rhoHat;
    initialRho(activeIndex) = initial.rhoHat;
    initialProjectionScore(activeIndex) = initial.projectionScore;
    finalProjectionScore(activeIndex) = final.projectionScore;
    initialResidualCost(activeIndex) = initial.residualCost;
    finalResidualCost(activeIndex) = final.residualCost;
    exitflag(activeIndex) = exitflagCurrent;
    iterations(activeIndex) = output.iterations;
    funcCount(activeIndex) = output.funcCount;
    optimizerObjective(activeIndex) = objectiveFinal;
    optimizerParameter(activeIndex, :) = qHat(:).';
end

assert(all(isfinite(delayFinal)));
assert(all(isfinite(xiFinal)));
assert(all(isfinite(zetaFinal)));
assert(all(isfinite(real(rhoFinal))));
assert(all(isfinite(imag(rhoFinal))));
assert(all(isfinite(initialProjectionScore)));
assert(all(isfinite(finalProjectionScore)));
assert(all(isfinite(initialResidualCost)));
assert(all(isfinite(finalResidualCost)));

result.activeRisIndices = activeRisIndices;
result.delayEstimate = delayFinal;
result.xiEstimate = xiFinal;
result.zetaEstimate = zetaFinal;
result.rhoEstimate = rhoFinal;
result.initialDelayEstimate = delayEstimate;
result.initialXiEstimate = xiEstimate;
result.initialZetaEstimate = zetaEstimate;
result.initialRhoEstimate = initialRho;
result.initialProjectionScore = initialProjectionScore;
result.finalProjectionScore = finalProjectionScore;
result.initialResidualCost = initialResidualCost;
result.finalResidualCost = finalResidualCost;
result.exitflag = exitflag;
result.iterations = iterations;
result.funcCount = funcCount;
result.optimizerObjective = optimizerObjective;
result.optimizerParameter = optimizerParameter;
result.options = options;
end

function objectiveValue = negativeProjectionScore( ...
    q, y, yEnergy, W, subcarrierIndex, cfg)
% q(3) 是无量纲 delay nu；优化过程中不做周期 wrap。
delayVector = exp(-1i * 2 * pi * subcarrierIndex * q(3));
a = upa_steering(q(1), q(2), cfg);
s = W.' * a;
M = delayVector * s.';
mu = M(:);

numerator = abs(mu' * y)^2;
muEnergy = real(mu' * mu);
projectionScore = numerator / (muEnergy * yEnergy);
objectiveValue = -projectionScore;
end

function evaluation = evaluatePhysicalParameters( ...
    y, W, xi, zeta, tau, subcarrierIndex, cfg)
nu = cfg.ofdm.subcarrier_spacing * tau;
delayVector = exp(-1i * 2 * pi * subcarrierIndex * nu);
a = upa_steering(xi, zeta, cfg);
s = W.' * a;
M = delayVector * s.';
mu = M(:);

muEnergy = real(mu' * mu);
yEnergy = real(y' * y);
rhoHat = (mu' * y) / (mu' * mu);
residual = y - rhoHat * mu;
residualCost = real(residual' * residual);
projectionScore = abs(mu' * y)^2 / (muEnergy * yEnergy);

evaluation.rhoHat = rhoHat;
evaluation.residualCost = residualCost;
evaluation.projectionScore = projectionScore;
end
