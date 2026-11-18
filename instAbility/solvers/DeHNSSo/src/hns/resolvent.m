function Res = resolvent(StabGrid, Opt, LHS, spec)
%% License (GNU GENERAL PUBLIC LICENSE v3)
%                  Delft Harmonic Navier-Stokes Solver
%     Copyright (C) 2026 P.C. Boldini, Y. Chen & M. Kotsonis
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
%     Please cite this code and the paper if you have used this for your
%     publication:

% DOI:

%% Description
% Resolvent analysis: the body force that gives the largest response.
%
%   Res = resolvent(StabGrid, Opt, LHS, spec)
%
%   A body force f on the momentum equations of mode (omega, beta) gives
%   ML q = B f. The gain is the energy of the response over the energy of the
%   forcing, sigma^2 = (q^H C^H W C q) / (f^H W f), with B and C the forcing
%   and response windows and W the quadrature weights. The optimal forcings
%   are the singular vectors of R = W_r^{1/2} C ML^{-1} B W_f^{-1/2}, computed
%   with svds from the LU factors of the linear run. Incompressible only.
%   Reference: Sipp & Marquet, Theor. Comput. Fluid Dyn. 27, 617 (2013),
%   doi:10.1007/s00162-012-0265-y.
%
%   spec fields:
%     omega, beta      mode (beta default 0)
%     k                number of modes (default 2)
%     xf, xr, yf, yr   forcing and response windows (default: up to the buffer)
%     comp_f, comp_r   components forced and measured (default [1 2] = u, v)
%     tol, maxit       svds tolerance and iteration cap (default 1e-8, 300)
%
%   Example, after a linear run with Opt.lu_mode = 'full':
%     [~, StabGrid, LHS] = main(Stab, StabGrid, Opt);
%     Res = resolvent(StabGrid, Opt, LHS, struct('omega', Stab.omega_0, 'k', 2));
%
%   Drop Res.F (the LU factors) before saving Res.

%% Inputs
% Name     size         unit   explanation
% StabGrid (1,1) struct [-]    stability grid, coordinate metrics and base-flow fields on it
% Opt      (1,1) struct [-]    solver options of the linear run
% LHS      (1,1) struct [-]    operator blocks and LU cache returned by main
% spec     (1,1) struct [-]    mode, number of modes and windows

%% Outputs
% Name size         unit   explanation
% Res  (1,1) struct [-]    gains sigma (1,k); forcing modes uf/vf/wf and response modes ur/vr/wr/pr
%                          (ny,nx,k), forcing of unit energy, response of energy sigma^2; the same as
%                          state vectors f/q; weights w; rows idx_f/idx_r; windows; svds flag; Vf, Rfun, F

ny = StabGrid.ny;  nx = StabGrid.nx;  nvar = 4;  N = nvar*ny*nx;
assert(size(LHS.MLA1,1) == N, 'resolvent: incompressible operator expected.');
assert(isfield(spec,'omega'), 'resolvent: spec.omega is required.');

omega  = spec.omega;
beta   = getfield_default(spec, 'beta',   0);
k      = getfield_default(spec, 'k',      2);
comp_f = getfield_default(spec, 'comp_f', [1 2]);
comp_r = getfield_default(spec, 'comp_r', [1 2]);
tol    = getfield_default(spec, 'tol',    1e-8);
maxit  = getfield_default(spec, 'maxit',  300);

%% Windows
x = StabGrid.x(:).';  y = StabGrid.y(:);
i_last = nx;
if isfield(Opt,'buffer') && any(strcmpi(Opt.buffer, {'on','para'}))
    i_last = round(Opt.xb*nx/100) - 1;
end
xf = getfield_default(spec, 'xf', [x(2) x(i_last)]);
xr = getfield_default(spec, 'xr', [x(2) x(i_last)]);
yf = getfield_default(spec, 'yf', [min(y) max(y)]);
yr = getfield_default(spec, 'yr', [min(y) max(y)]);

idx_f = window_rows(x, y, xf, yf, comp_f, i_last, nvar);
idx_r = window_rows(x, y, xr, yr, comp_r, i_last, nvar);

%% Quadrature weights (area per node)
J = ones(ny, nx);
if isfield(StabGrid, 'xix')
    J = abs(StabGrid.xix.*StabGrid.etay - StabGrid.xiy.*StabGrid.etax);
end
w  = (trapz_weights(y) * trapz_weights(x).') ./ J;
wN = repmat(reshape(w, ny, 1, nx), [1 nvar 1]);
wN = wN(:);
sf = sqrt(wN(idx_f));  sr = sqrt(wN(idx_r));

%% LU of the operator: from the linear run, else factorise as do_lu
key = matlab.lang.makeValidName(sprintf('o%g_b%g', omega, beta));
if isfield(LHS,'LU_cache') && isfield(LHS.LU_cache, key)
    F = LHS.LU_cache.(key);
else
    ML = assemble_ML(omega, beta, LHS.MLA1, LHS.MLA2, LHS.MLA3, LHS.MLA4, LHS.MLA5, LHS.MFD_diff);
    fprintf(' Resolvent: factorising %s (%d unknowns) ...', key, N);  t0 = tic;
    if isfield(Opt,'lu_ordering') && strcmpi(Opt.lu_ordering, 'symrcm')
        r = symrcm(ML);
        [L, U, p] = lu(ML(r,r), 'vector');
        F = struct('L', L, 'U', U, 'p', p, 'r', r);
    else
        [L, U, p, q] = lu(ML, 'vector');
        F = struct('L', L, 'U', U, 'p', p, 'q', q);
    end
    fprintf(' %.1f s\n', toc(t0));
end

%% Singular values
Rfun = @(v, flag) apply_R(v, flag, F, idx_f, idx_r, sf, sr, N);
fprintf(' Resolvent: svds for %d mode(s), %d forced and %d measured rows ...', k, numel(idx_f), numel(idx_r));
t0 = tic;
[~, S, Vf, flag] = svds(Rfun, [numel(idx_r) numel(idx_f)], k, 'largest', ...
                        struct('Tolerance', tol, 'MaxIterations', maxit));
sigma = diag(S).';
fprintf(' %.1f s (flag %d)\n', toc(t0), flag);

%% Modes on the grid
f = zeros(N, k);  q = zeros(N, k);
for m = 1:k
    f(idx_f, m) = Vf(:, m) ./ sf;
    q(:, m)     = lu_solve(F, f(:, m), false);
    [~, i] = max(abs(q(:, m)));                       % phase: largest response entry real
    ph = exp(-1i*angle(q(i, m)));
    f(:, m) = f(:, m)*ph;  q(:, m) = q(:, m)*ph;
end
Fm = reshape(f, ny, nvar, nx, k);  Qm = reshape(q, ny, nvar, nx, k);

Res = struct('sigma', sigma, ...
             'uf', squeeze(Fm(:,1,:,:)), 'vf', squeeze(Fm(:,2,:,:)), 'wf', squeeze(Fm(:,3,:,:)), ...
             'ur', squeeze(Qm(:,1,:,:)), 'vr', squeeze(Qm(:,2,:,:)), 'wr', squeeze(Qm(:,3,:,:)), ...
             'pr', squeeze(Qm(:,4,:,:)), 'f', f, 'q', q, 'w', w, 'idx_f', idx_f, 'idx_r', idx_r, ...
             'xf', xf, 'xr', xr, 'yf', yf, 'yr', yr, 'omega', omega, 'beta', beta, ...
             'flag', flag, 'Vf', Vf, 'Rfun', Rfun, 'F', F);
fprintf(' Resolvent gains: %s\n', mat2str(sigma, 6));
end

%% ---------------------------------------------------------------------------
function v = getfield_default(s, name, default)
v = default;
if isfield(s, name) && ~isempty(s.(name)); v = s.(name); end
end

function idx = window_rows(x, y, xw, yw, comps, i_last, nvar)
% Rows of the given components inside the window, without station 1, the
% buffer and the boundary rows.
ny = numel(y);  nx = numel(x);
in_x = false(1, nx);  in_x(2:i_last) = x(2:i_last) >= xw(1) & x(2:i_last) <= xw(2);
in_y = y >= yw(1) & y <= yw(2);  in_y([1 ny]) = false;
mask = false(ny, nvar, nx);
mask(in_y, comps, in_x) = true;
idx  = find(mask(:));
assert(~isempty(idx), 'resolvent: empty window.');
end

function w = trapz_weights(s)
d = abs(diff(s(:)));
w = [d; 0]/2 + [0; d]/2;
end

function y = apply_R(v, flag, F, idx_f, idx_r, sf, sr, N)
% R v and R^H v.
if strcmp(flag, 'notransp')
    b = zeros(N, 1);  b(idx_f) = v./sf;
    q = lu_solve(F, b, false);
    y = sr.*q(idx_r);
else
    b = zeros(N, 1);  b(idx_r) = sr.*v;
    g = lu_solve(F, b, true);
    y = g(idx_f)./sf;
end
end

function x = lu_solve(F, b, transp)
% ML x = b, or ML^H x = b if transp, for the factor forms of do_lu.
if isfield(F, 'r');  b = b(F.r);  end
x = zeros(size(b));
if isfield(F, 'D')
    if ~transp;  x = F.Q*(F.U\(F.L\(F.P*(F.D\b))));
    else;        x = F.D'\(F.P'*(F.L'\(F.U'\(F.Q'*b))));  end
elseif isfield(F, 'q')
    if ~transp;  x(F.q) = F.U\(F.L\b(F.p));
    else;        x(F.p) = F.L'\(F.U'\b(F.q));  end
elseif isfield(F, 'L')
    if ~transp;  x = F.U\(F.L\b(F.p));
    else;        x(F.p) = F.L'\(F.U'\b);  end
elseif isfield(F, 'dML')
    if ~transp;  x = F.dML\b;  else;  x = F.dML'\b;  end
else
    if ~transp;  x = F\b;  else;  x = F'\b;  end
end
if isfield(F, 'r');  x(F.r) = x;  end
end
