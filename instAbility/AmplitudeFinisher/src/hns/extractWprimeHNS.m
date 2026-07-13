function npse = extractWprimeHNS(out, in)
% EXTRACTWPRIMEHNS  DeHNSSo StabRes -> w' mode RMS profiles + combined RMS.
%
%   npse = extractWprimeHNS(out, in)
%
% From a runHNS bundle, extract the first nModes spanwise Fourier modes of w'
% and reduce them to the combined RMS profile that is the match target.
%
% After postprocess, StabRes stores (conjugates removed):
%   betavec  1 x nmode         wavenumbers [0, beta_0, 2*beta_0, ...]
%   u/v/w/p  nmode x ny x nx    NORMALIZED shape functions (|u_max| = 1)
%   A        nmode x nx         per-mode amplitude (A = 2|u_max| for non-MFD)
% so the physical spanwise-Fourier mode field is w(idx)*A(idx)/2, and its
% spanwise RMS carries a convention factor in.match.rmsFactor.
%
% Modes are selected by physical wavenumber (betavec == k*beta_0), robust to
% ordering. Coordinates are re-dimensionalized by lref to metres, arc-length
% from the LE (+ in.match.xOffset to align with the PIV frame).
%
% Returns:
%   npse.x     1 x nx        streamwise stations [m]
%   npse.y     ny x 1        wall-normal coordinate [m]
%   npse.mode  ny x nx x K   per-mode w' RMS profiles [m/s]
%   npse.rms   ny x nx       combined RMS sqrt(sum_k mode_k^2) [m/s]

res  = out.res;
g    = out.grid;
lref = out.lref;
Uref = out.Uref;                 % [m/s] StabRes.w is non-dim (velocity/Uref)
b0   = 2*pi * lref / in.Stab.lambda_z;
tol  = 1e-6 * max(1, b0);

if ~isfield(res,'betavec') || ~isfield(res,'w')
    error('extractWprimeHNS:badOutput','StabRes lacks betavec/w.');
end

% Number of spanwise harmonics to extract for the diagnostics: matches the
% PIV modes present (set from loadPIVwPrime), capped by what the HNS solve
% actually contains; falls back to the available harmonics if unset.
nKmax = sum(res.betavec > tol);
if isfield(in.validation,'nModes') && ~isempty(in.validation.nModes)
    nK = min(in.validation.nModes, nKmax);
else
    nK = nKmax;
end

% --- coordinates (physical, from LE) ------------------------------------
npse.y = local_phys(g.etaun(:),  lref);              % ny x 1  [m]
npse.x = local_phys(g.x(:).',    lref) + in.match.xOffset;   % 1 x nx  [m]

[nmode, ny, nx] = size(res.w);                       %#ok<ASGLU>

% Number of spanwise harmonics present (beta = k*b0, k >= 1; excludes the MFD
% k=0, which is the spanwise-mean distortion and NOT part of w').
Nsp = sum(res.betavec > tol);

phys = @(idx) squeeze(abs(res.w(idx,:,:))) .* (res.A(idx,:)/2) ...
              .* in.match.rmsFactor .* Uref;          % ny x nx  [m/s]

% --- per-mode w' RMS profiles (first nK, for diagnostics) ----------------
npse.mode = zeros(ny, nx, nK);
for k = 1:nK
    idx = find(abs(res.betavec - k*b0) <= tol, 1, 'first');
    if isempty(idx)
        error('extractWprimeHNS:noMode', ...
              'No HNS mode with beta = %.4g (k=%d). betavec=[%s]', ...
              k*b0, k, num2str(res.betavec));
    end
    npse.mode(:,:,k) = phys(idx);
end
npse.rms = sqrt(sum(npse.mode.^2, 3));                % per-mode RMS (diagnostic, first nK)

% --- full-plane RMS over ALL spanwise harmonics (peak-match target) ------
% By Parseval, the spanwise RMS of the reconstructed y-z plane is
% sqrt(sum_{k=1}^{Nsp} |w'_k|^2) — matches the raw PIV plane RMS definition.
sumsq = zeros(ny, nx);
for k = 1:Nsp
    idx = find(abs(res.betavec - k*b0) <= tol, 1, 'first');
    if ~isempty(idx); sumsq = sumsq + phys(idx).^2; end
end
npse.rmsFull = sqrt(sumsq);
end

% -------------------------------------------------------------------------
function p = local_phys(c, lref)
% Return the coordinate in metres. StabGrid coordinates are non-dimensional
% (by lref) when their magnitude exceeds 1; if already O(<1) they are metres.
if max(abs(c(:))) > 1
    p = c * lref;
else
    p = c;
end
end
