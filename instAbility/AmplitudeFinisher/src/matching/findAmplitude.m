function result = findAmplitude(StabGrid, piv, in, seed)
% FINDAMPLITUDE  Pin the disturbance amplitude Stab.A0_fund to the PIV.
%
%   result = findAmplitude(StabGrid, piv, in)
%
% Strategy (the "finisher"):
%   1. LINEAR seed  — one linear HNS solve excites only the fundamental; its
%      w'-RMS profile is projected onto the PIV target to get a seed amplitude.
%   2. NONLINEAR bracket sweep — solve at A0 = seed .* in.match.sweep (fixed N),
%      score each against the PIV, and interpolate the A0 where the mean HNS
%      peak equals the PIV peak (peak-ratio -> 1) — the physical CFI-peak
%      amplitude.
%
% Returns:
%   result.A0            best amplitude (peak-interpolated for peakRMS)
%   result.A0sweptBest   best amplitude actually solved in the sweep
%   result.J             objective at result.A0sweptBest
%   result.N             mode count
%   result.hns           extracted HNS w' at result.A0sweptBest
%   result.seed          linear-seed amplitude
%   result.sweep         per-point struct (A0, J, ok, ratio, hns)

N = in.Stab.N;

%% 1. Amplitude seed
% Prefer the linear-check A0_fund (anchored at in.linearCheck.xcAnchor) passed in
% by the caller -- one linear estimate feeds BOTH the linear-check figure and the
% centre of this sweep. Fall back to a window-projected linear seed if not given.
if nargin >= 4 && ~isempty(seed) && isfinite(seed) && seed > 0
    fprintf('\n[seed] linear-check A0_fund = %.4g (anchored x/c = %g%%)\n', ...
            seed, in.linearCheck.xcAnchor);
else
    fprintf('\n[seed] linear solve ...\n');
    outL = runHNS(StabGrid, in.Stab.A0_fund, N, in, 'linear');
    hnsL = extractWprimeHNS(outL, in);
    seed = local_linear_seed(hnsL, outL, piv, in);
    fprintf('[seed] A0_seed = %.4g\n', seed);
end

%% 2. Nonlinear bracket sweep
A0list = seed * in.match.sweep;
Sw = struct('A0',{},'J',{},'ok',{},'ratio',{},'hns',{});
for q = 1:numel(A0list)
    A0 = A0list(q);
    fprintf('\n[sweep %d/%d] A0 = %.4g (N=%d) ...\n', q, numel(A0list), A0, N);
    t = tic;
    try
        out = runHNS(StabGrid, A0, N, in, 'nonlinear');
        hns = extractWprimeHNS(out, in);
        ok  = all(isfinite(out.res.A(:)));
        J   = wPrimeObjective(hns, piv, in);
        r   = local_peak_ratio(hns, piv, in);
        Sw(end+1) = struct('A0',A0,'J',J,'ok',ok,'ratio',r,'hns',hns); %#ok<AGROW>
        fprintf('   %.0fs | ok=%d | peak HNS/PIV = %.3f | J = %.4g\n', toc(t), ok, r, J);
    catch ME
        fprintf('   FAILED: %s\n', ME.message);
    end
end

%% 3. Pick the amplitude
good = [Sw.ok] & isfinite([Sw.J]);
if ~any(good)
    error('findAmplitude:noConverged', ...
          'No converged sweep point — try a lower amplitude window or larger N.');
end
Sg = Sw(good);
[~, ib] = min([Sg.J]);   best = Sg(ib);

A0 = best.A0;
if numel(Sg) >= 2
    % Convention: return the A0 at which the MEAN peak ratio over the window,
    % mean_i( HNS_peak_i / PIV_peak_i ), equals 1 (interpolated from the sweep).
    % This centers the per-station ratios about 1 (over/under-shoot cancel) — it
    % is NOT the minimizer of sum_i |ratio_i - 1| (that is A0sweptBest/J below).
    A0 = interp1([Sg.ratio], [Sg.A0], 1.0, 'linear', 'extrap');
end

result = struct('A0', A0, 'A0sweptBest', best.A0, 'J', best.J, 'N', N, ...
                'hns', best.hns, 'seed', seed, 'sweep', Sw);
end

% ---------------------------------------------------------------------
function seed = local_linear_seed(hnsL, outL, piv, in)
% Project the linear w'-RMS profile onto the PIV target over the scored
% stations, then map the fit scale to an inflow fundamental amplitude.
sel = local_scored(piv, in);
num = 0; den = 0;
for i = sel
    [~, ix] = min(abs(hnsL.x - piv.x(i)));
    f = interp1(hnsL.y, hnsL.rmsFull(:, ix), piv.y{i}, 'linear', 'extrap');
    g = local_target(piv, i);
    num = num + sum(f .* g);   den = den + sum(f .* f);
end
c = num / max(den, eps);
b0 = 2*pi * outL.lref / in.Stab.lambda_z;
fidx = find(abs(outL.res.betavec - b0) <= 1e-6*max(1,b0), 1, 'first');
seed = c * outL.res.A(fidx, 1) / 2;
seed = min(max(seed, 1e-6), 1e-1);
end

% ---------------------------------------------------------------------
function r = local_peak_ratio(hns, piv, in)
% Mean over scored stations of (HNS main-peak / PIV main-peak).
sel  = local_scored(piv, in);
rr = zeros(1, numel(sel));
for k = 1:numel(sel)
    i = sel(k);
    [~, ix] = min(abs(hns.x - piv.x(i)));
    yq = piv.y{i};
    hn = interp1(hns.y, hns.rmsFull(:, ix), yq, 'linear', 'extrap');
    tp = local_target(piv, i);
    rr(k) = max(hn) / max(tp);
end
r = mean(rr);
end

% ---------------------------------------------------------------------
function g = local_target(piv, i)
% PIV target profile: raw plane RMS if available, else the mode-sum.
if isfield(piv, 'rmsFull') && ~isempty(piv.rmsFull{i})
    g = piv.rmsFull{i};
else
    g = piv.rms{i};
end
end

% ---------------------------------------------------------------------
function sel = local_scored(piv, in)
if ~isempty(in.match.stationIdx)
    sel = in.match.stationIdx(:).';
else
    sel = 1:numel(piv.x);
end
end
