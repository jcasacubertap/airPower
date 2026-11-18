function BF_new = resample_cartesian(BF, params)
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
% Interpolate a structured BF onto a Malik-Chebyshev η × equidistant ξ grid.
%
%   BF_new = resample_cartesian(BF, params)
%
%   For BF data already on a structured (nx_raw × ny_raw) grid but at a different
%   resolution than DeHNSSo expects; resample_curvilinear and
%   resample_unstructured handle the body-fitted equivalents.
%
%   BF holds X, Y (1-D axes), U, V, W on a regular grid, or a flat scatter
%   canonicalised on entry. params: n_eta_new / n_xi_new (default source);
%   H, y_i (0.99 span, H/10); xi_range / xi_trim_* to restrict or trim the span;
%   mug, sig, ag and spread for the streamwise clustering (see cluster_xi);
%   interp_method (default 'makima'); rescale, Uref, lref to non-dimensionalise;
%
%   Output: BF_new (X, Y, U, V, W, H, y_i, cluster) plus any Re/lref/Uref/nu.

%% Inputs
% Name   size         unit   explanation
% BF     (1,1) struct [-]    base flow on its own grid
% params (1,1) struct [-]    grid-generation parameters

%% Outputs
% Name   size         unit   explanation
% BF_new (1,1) struct [-]    base flow resampled onto the stability grid

%% Canonicalise a flat-scatter input to gridded (ny x nx); no-op if already 2-D
BF = reshape_scatter(BF);

x_raw  = BF.X(:).';          % [1 × nx_raw]
y_raw  = BF.Y(:);            % [ny_raw × 1]
nx_raw = numel(x_raw);
ny_raw = numel(y_raw);
U_raw  = BF.U;   V_raw = BF.V;   W_raw = BF.W;
% Compressible fields (auto-detected by presence of ρ and T).
is_compr = isfield(BF,'rho') && isfield(BF,'T');
if is_compr
    rho_raw = BF.rho;
    T_raw   = BF.T;
end

%% Optional non-dim rescale: U/Uref, X/lref, Y/lref
do_rescale = isfield(params, 'rescale') && params.rescale;
if do_rescale
    assert(isfield(params,'Uref') && isfield(params,'lref'), ...
        'resample_cartesian: rescale=true requires params.Uref and params.lref.');
    U_raw = U_raw / params.Uref;
    V_raw = V_raw / params.Uref;
    W_raw = W_raw / params.Uref;
    x_raw = x_raw / params.lref;
    y_raw = y_raw / params.lref;
    fprintf('Rescaled: U/=%.3e (Uref),  X,Y/=%.3e (lref)\n', params.Uref, params.lref);
end

%% Target grid dimensions — default to source resolution
if isfield(params, 'n_eta_new') && ~isempty(params.n_eta_new)
    ny_new = params.n_eta_new;
else
    ny_new = ny_raw;
end

if isfield(params, 'n_xi_new') && ~isempty(params.n_xi_new)
    nx_new = params.n_xi_new;
else
    nx_new = nx_raw;
end

%% Wall-normal domain (H, y_i)
if isfield(params, 'H') && ~isempty(params.H)
    H = params.H;
elseif isfield(BF, 'H') && ~isempty(BF.H)
    H = BF.H;
else
    H = abs(y_raw(end) - y_raw(1)) * 0.99;   % BF.Y may be descending (η: top -> wall)
end
if isfield(params, 'y_i') && ~isempty(params.y_i)
    y_i = params.y_i;
elseif isfield(BF, 'y_i') && ~isempty(BF.y_i)
    y_i = BF.y_i;
else
    y_i = H / 10;
end
if y_i >= H/2
    y_i = H/3;
    warning('y_i clamped to H/3 = %.5f (must be < H/2 for Malik mapping).', y_i);
end

%% Malik-Chebyshev η — descending: η(1)=H (freestream), η(end)=0 (wall)
a = H * y_i / (H - 2*y_i);
b = 1 + 2*a/H;
k_vec   = (0 : ny_new - 1)';
z_cheb  = cos(pi * k_vec / (ny_new - 1));
eta_new = a * (1 + z_cheb) ./ (b - z_cheb);             % [ny_new × 1], η in [0, H]

%% Streamwise ξ — optionally trimmed, optionally clustered
% Priority: xi_range overrides xi_trim_*. xi_trim_* drops a fraction from
% each end of the x-span (mirrors curvilinear convention).
if isfield(params, 'xi_range') && ~isempty(params.xi_range)
    x_lo = x_raw(1) + params.xi_range(1) * (x_raw(end) - x_raw(1));
    x_hi = x_raw(1) + params.xi_range(2) * (x_raw(end) - x_raw(1));
else
    trim_lo = 0;  trim_hi = 0;
    if isfield(params, 'xi_trim_inflow')  && ~isempty(params.xi_trim_inflow)
        trim_lo = params.xi_trim_inflow;
    end
    if isfield(params, 'xi_trim_outflow') && ~isempty(params.xi_trim_outflow)
        trim_hi = params.xi_trim_outflow;
    end
    x_lo = x_raw(1) + trim_lo       * (x_raw(end) - x_raw(1));
    x_hi = x_raw(end) - trim_hi     * (x_raw(end) - x_raw(1));
end

% xi is x here, so the range ends map directly onto the clustering axis.
% spread (refinement half-width) sizes sig to cover ~3x that width.
mug   = get_opt(params, 'mug', []);
ag    = get_opt(params, 'ag',  1);
spread = get_opt(params, 'spread', []);
if ~isempty(spread)
    sig = 3 * spread ./ (0.5 * (x_hi - x_lo));
else
    sig = get_opt(params, 'sig', 0.2);
end
xi_new = cluster_xi(x_lo, x_hi, nx_new, mug, sig, ag, x_lo, x_hi);   % [1 × nx_new]
report_clustering(mug, sig, ag, xi_new);

%% Solid: put a station on the face and a node on the top
% Both are uniform rescales, as in DeHNSSo-v1: the Malik map scales exactly with
% (H, y_i), so the node distribution keeps its shape and the wall location is
% exact without ghost rows.
step = get_opt(params, 'step', struct('type','FLAT'));
if ~isfield(step, 'seg'); step.seg = zeros(0,3); end
step.snapped = false;
% Two rescales place one face and one height, so only a single box can be put on
% grid lines. 
snap = get_opt(params, 'snap', true);
if size(step.seg,1) > 1
    fprintf('Solid: %d boxes, so the grid is not snapped and the ghost-cell boundary carries the surfaces\n', size(step.seg,1));
elseif ~snap && ~strcmpi(step.type, 'FLAT')
    fprintf('Solid: not snapped (params.snap = false), the ghost-cell boundary carries the surface\n');
elseif ~strcmpi(step.type, 'FLAT') && step.x > x_lo && step.x < x_hi && step.h > 0 && step.h < H
    [~, k] = min(abs(xi_new - step.x));
    xi_new = x_lo + (xi_new - x_lo) * (step.x - x_lo) / (xi_new(k) - x_lo);
    [~, j] = min(abs(eta_new - step.h));
    f_eta   = step.h / eta_new(j);
    eta_new = f_eta * eta_new;  H = f_eta * H;  y_i = f_eta * y_i;
    step.i = k;  step.j = j;  step.snapped = true;
    fprintf('Solid: station %d snapped to the face x = %.4f (domain end %.4f -> %.4f), node %d to the top h = %.4f (H scaled by %.5f)\n', ...
            k, step.x, x_hi, xi_new(end), j, step.h, f_eta);
end

%% Interpolate onto the target grid ('makima': no overshoot on non-uniform data)
interp_method = 'makima';
if isfield(params, 'interp_method') && ~isempty(params.interp_method)
    interp_method = params.interp_method;
end
[Xq, Yq] = meshgrid(xi_new, eta_new);
U_new = interp2(x_raw, y_raw, U_raw, Xq, Yq, interp_method);
V_new = interp2(x_raw, y_raw, V_raw, Xq, Yq, interp_method);
W_new = interp2(x_raw, y_raw, W_raw, Xq, Yq, interp_method);
% ρ, T (already non-dim by the BL caller; not affected by params.rescale)
if is_compr
    rho_new = interp2(x_raw, y_raw, rho_raw, Xq, Yq, interp_method);
    T_new   = interp2(x_raw, y_raw, T_raw,   Xq, Yq, interp_method);
end

%% Package into DeHNSSo-ready BF (no mesh_type -> build_stab_grid Cartesian branch)
BF_new       = struct();
BF_new.X     = xi_new;
BF_new.Y     = eta_new;
BF_new.U     = U_new;
BF_new.V     = V_new;
BF_new.W     = W_new;
BF_new.H     = H;
BF_new.y_i   = y_i;
BF_new.xi1D  = xi_new;
BF_new.eta1D = eta_new;

% Clustering metadata for the caller (xi = x, so xc_xi = mug); w is the
% provisioned actuator half-width. See build_stab_grid / add_wall_forcing.
BF_new.cluster = struct('mug', mug, 'sig', sig, 'ag', ag, 'xc_xi', mug, 'w', spread);

% Compressible fields + thermo metadata (build_stab_grid runs Sutherland).
if is_compr
    BF_new.rho = rho_new;
    BF_new.T   = T_new;
    for f_name = {'Ma','Pr','gamma','Sstar'}
        if isfield(BF, f_name{1}) && ~isempty(BF.(f_name{1}))
            BF_new.(f_name{1}) = BF.(f_name{1});
        end
    end
end

for dn = bf_deriv_fields()
    if isfield(BF, dn{1}) && ~isempty(BF.(dn{1}))
        BF_new.(dn{1}) = interp2(x_raw, y_raw, BF.(dn{1}), Xq, Yq, interp_method);
    end
end

% Interpolation smears the solid edges, so restore the zeros inside every box.
if ~strcmpi(step.type, 'FLAT')
    in_solid = false(size(Xq));
    for b = 1:size(step.seg, 1)
        xb = step.seg(b,1);  xe = step.seg(b,2);  hb = step.seg(b,3);
        in_solid = in_solid | (Xq >= xb - 1e-9*max(1,abs(xb)) & Xq <= xe + 1e-9*max(1,abs(xe)) ...
                             & Yq > 0 & Yq <= hb + 1e-9*hb);
    end
    for f_name = [{'U','V','W'}, bf_deriv_fields()]
        if isfield(BF_new, f_name{1});  BF_new.(f_name{1})(in_solid) = 0;  end
    end
end
BF_new.step = step;

fprintf('Cartesian resample done: %d × %d → %d × %d (Malik-Cheb η, %s interp)\n', ...
        ny_raw, nx_raw, ny_new, nx_new, interp_method);
end


function BF = reshape_scatter(BF)
%RESHAPE_SCATTER  Flat-lattice scatter -> canonical (ny × nx) gridded form.
if ~isvector(BF.X) || ~isvector(BF.Y) || numel(BF.X) ~= numel(BF.Y)
    return;
end

tolX = max(1e-6 * (max(BF.X) - min(BF.X)), 1e-12);
tolY = max(1e-6 * (max(BF.Y) - min(BF.Y)), 1e-12);
[Xu, ~, ix] = uniquetol(BF.X(:), tolX, 'DataScale', 1);
[Yu, ~, iy] = uniquetol(BF.Y(:), tolY, 'DataScale', 1);
nx = numel(Xu);
ny = numel(Yu);

fprintf('Reshaping flat scatter → (ny=%d × nx=%d)\n', ny, nx);

fnames = setdiff(fieldnames(BF), {'X','Y'});
for k = 1:numel(fnames)
    f   = fnames{k};
    val = BF.(f);
    if ~isnumeric(val) || numel(val) ~= numel(ix)
        continue
    end
    M   = nan(ny, nx);
    lin = sub2ind([ny, nx], iy, ix);
    M(lin) = val(:);
    if any(isnan(M(:)))
        warning('resample_cartesian:gaps', ...
            'Field %s has %d missing lattice points after reshape.', f, nnz(isnan(M)));
    end
    BF.(f) = M;
end

BF.X = Xu(:).';   % 1 × nx
BF.Y = Yu(:);     % ny × 1
end
