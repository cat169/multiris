function [result, model] = compute_channel_fim( ...
    qHat, rhoHat, W, separatedNoiseVariance, cfg)
%COMPUTE_CHANNEL_FIM Equivalent information after profiling complex gain.
% qHat = [xi; zeta; nu], nu = cfg.ofdm.subcarrier_spacing * tau.
% Noise convention: E[abs(n).^2] = separatedNoiseVariance AFTER DFT
% separation. The factor 2 below is for real parameters in complex noise.
% W is N-by-G, with the same nonconjugate transpose as channel MLE.
% Optional second output exposes the atom and analytic derivatives for
% lightweight validation. No inverse, clipping or regularization is used.

assert(isnumeric(qHat) && isreal(qHat) && ...
    isequal(size(qHat), [3, 1]) && all(isfinite(qHat)), ...
    'compute_channel_fim:InvalidParameters', ...
    'qHat must be a finite real 3-by-1 vector [xi; zeta; nu].');
assert(isnumeric(rhoHat) && isscalar(rhoHat) && isfinite(rhoHat), ...
    'compute_channel_fim:InvalidGain', 'rhoHat must be a finite scalar.');
assert(isnumeric(separatedNoiseVariance) && ...
    isreal(separatedNoiseVariance) && isscalar(separatedNoiseVariance) && ...
    isfinite(separatedNoiseVariance) && separatedNoiseVariance > 0, ...
    'compute_channel_fim:PositiveNoiseRequired', ...
    'Fisher information requires positive separated complex noise variance.');
assert(isnumeric(W) && ismatrix(W) && ~isempty(W) && ...
    size(W, 1) == cfg.ris.nx * cfg.ris.ny && all(isfinite(W(:))), ...
    'compute_channel_fim:InvalidReflections', ...
    'W must be finite and have nx*ny rows.');

indexX = (0:(cfg.ris.nx - 1)).';
indexY = (0:(cfg.ris.ny - 1)).';
phaseX = (2 * pi / cfg.constants.lambda) * cfg.ris.dx * indexX;
phaseY = (2 * pi / cfg.constants.lambda) * cfg.ris.dy * indexY;
ax = exp(1i * phaseX * qHat(1));
ay = exp(1i * phaseY * qHat(2));
a = kron(ax, ay);                  % x slow, y fast: upa_steering convention
aXi = kron(1i * phaseX .* ax, ay);
aZeta = kron(ax, 1i * phaseY .* ay);
k = (0:(cfg.ofdm.num_subcarriers - 1)).';
d = exp(-1i * 2 * pi * k * qHat(3));
dNu = -1i * 2 * pi * k .* d;
s = W.' * a;                      % deliberately NOT W' * a
mu = kron(s, d);                  % vec(d * s.')
D = [kron(W.' * aXi, d), kron(W.' * aZeta, d), kron(s, dNu)];
muEnergy = real(mu' * mu);
assert(isfinite(muEnergy) && muEnergy > 0, ...
    'compute_channel_fim:ZeroAtom', ...
    'The atom must have positive finite energy to eliminate complex gain.');

% Schur complement for nuisance [real(rho), imag(rho)]; no large projector.
crossTerm = D' * mu;
gram = D' * D - (crossTerm * crossTerm') / muEnergy;
fim = (2 * abs(rhoHat)^2 / separatedNoiseVariance) * real(gram);
fim = (fim + fim.') / 2;
assert(all(isfinite(fim(:))), 'compute_channel_fim:NonfiniteInformation', ...
    'Fisher information overflowed or contains nonfinite values.');

result.fim = fim;
result.eigenvalues = sort(eig(fim));
result.conditionNumber = cond(fim);
result.parameterOrder = {'xi'; 'zeta'; 'nu'};
% Relative roundoff allowance only; it never changes the matrix/eigenvalues.
result.psdTolerance = 100 * eps * max(norm(fim, 2), realmin);
result.isPositiveSemidefinite = ...
    result.eigenvalues(1) >= -result.psdTolerance;
result.isNearlySingular = rcond(fim) < sqrt(eps);
if nargout > 1
    model.atom = mu;
    model.derivatives = D;
end
end
