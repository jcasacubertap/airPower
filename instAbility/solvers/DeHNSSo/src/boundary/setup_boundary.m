function [Bound]=setup_boundary(Grid,Opt,Stab)
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
% Set up everything the operator needs to know about the boundaries.
%
%   Bound = setup_boundary(Grid, Opt, Stab)
%
%   Four independent pieces, each in its own block below:
%     buffer      - the outflow zone where the solution is damped ('on') or
%                   parabolised ('para'). The tanh ramp starts at Opt.xb percent
%                   of the domain and steepens with Opt.kappa; the nonlinear
%                   forcing gets its own ramp from Opt.nltbufxb.
%     wall masks  - a solid standing on the wall (Grid.StepType ~= 'FLAT'): the
%                   nodes inside it, the height profile, and either a ghost row
%                   per station or plain masking (Opt.solid).
%     inhomog. BC - Stab.bcw / Stab.bct interpolated from wherever the caller
%                   defined them onto the numerical grid.
%     resonator   - Opt.hr resolved to an admittance and a mouth mask.
%
%   Bound fields, all read downstream by LHS / RHS / solve_modes:
%     bufc, bufp, bufcf, ibuf, ibuff       buffer ramps and start indices
%     bufs, bufsp, h_surf                  solid masks (velocity, continuity), height
%     ghost_idx, ghost_wts, stencil_wall   ghost-cell row per station and its weights
%     bcw, bct                             (n_bc x nx x nf), n_bc = 3 | 4 (compressible)
%     hr                                   only when Opt.hr is set

%% Inputs
% Name size         unit   explanation
% Grid (1,1) struct [-]    grid sizes, coordinates and inverse metric terms
% Opt  (1,1) struct [-]    solver options; see the option index in docs/DeHNSSo_guide.md
% Stab (1,1) struct [-]    disturbance specification: modal truncation, frequencies, inflow condition and amplitudes

%% Outputs
% Name  size         unit   explanation
% Bound (1,1) struct [-]    boundary conditions and outflow-buffer coefficients

%% Buffer

switch lower(Opt.buffer)
case 'on'
    % Amplitude-damping buffer: scales the whole finite-difference (FD) block by bufc -> 0
    Bound.ibuf  = round(Opt.xb * Grid.nx / 100);
    Bound.bufc  = ones(1, Grid.nx);
    Bound.bufp  = ones(1, Grid.nx);   % Fi (streamwise diffusion term) is left unchanged in 'on' mode
    Bound.bufc  = buffer_ramp(Bound.bufc,  Bound.ibuf,  Grid.nx, Opt.kappa);

    % NLT forcing buffer
    Bound.ibuff = round(Opt.nltbufxb * Grid.nx / 100);
    Bound.bufcf = ones(1, Grid.nx);
    Bound.bufcf = buffer_ramp(Bound.bufcf, Bound.ibuff, Grid.nx, Opt.kappa);

    fprintf(' Buffer: amplitude damping from %.3g%% to 100%%\n', Opt.xb);

case 'para'
    % Parabolization buffer: ramps Fi (streamwise diffusion) to zero in the buffer,
    % keeping advection and wall-normal diffusion. This makes the outflow behave
    % like a linear parabolized stability equation (LPSE) solve, avoiding reflections.
    Bound.ibuf  = round(Opt.xb * Grid.nx / 100);
    Bound.bufc  = ones(1, Grid.nx);   % no amplitude damping in this mode

    Bound.bufp  = ones(1, Grid.nx);
    Bound.bufp  = buffer_ramp(Bound.bufp,  Bound.ibuf,  Grid.nx, Opt.kappa);

    % NLT forcing buffer (same ramp as bufp)
    Bound.ibuff = round(Opt.nltbufxb * Grid.nx / 100);
    Bound.bufcf = ones(1, Grid.nx);
    Bound.bufcf = buffer_ramp(Bound.bufcf, Bound.ibuff, Grid.nx, Opt.kappa);

    fprintf(' Buffer: parabolisation from %.3g%% to 100%%\n', Opt.xb);

case 'off'
    Bound.ibuf  = Grid.nx + 1;
    Bound.ibuff = Grid.nx + 1;
    Bound.bufc  = ones(1, Grid.nx);
    Bound.bufp  = ones(1, Grid.nx);
    Bound.bufcf = ones(1, Grid.nx);
    fprintf(' Buffer: off\n');
end

%% Wall masks: a solid standing on the wall
% A step, a staircase, a backward-facing step or a block on the wall is described
% by Grid.step.seg: one box [x_lo, x_hi, h] per run of constant height. Nothing
% here changes on a flat wall; the masks stay all ones and the solver never
% notices them.
%
% Two masks come out of this block. bufs multiplies the velocity rows and bufsp
% the continuity rows in LHS: a zero switches the row off, and Gi/Pres then
% impose q = 0 and dp/dy = 0 there instead. Where the top of the solid does not
% fall on a grid node, Opt.solid decides what to do:
%   'mask'   zero the nodes inside the solid and let the first masked node carry
%            the no-slip condition. The wall sits up to one cell too low.
%   'ghost'  keep that first node as a ghost point whose value is interpolated
%            (cubic Lagrange) so that the velocity is zero at the exact height
%            (Tseng & Ferziger 2003, JCP 192). Exact, but it rings at 2*dx at a
%            step face.
%   'auto'   'mask' when gridgen snapped a node onto the surface, so the wall is
%            exact anyway; 'ghost' otherwise.
Grid = step_geometry(Grid);                % accepts Grid.step or the bare scalars
Bound.bufs      = ones(Grid.ny, Grid.nx);  % velocity rows
Bound.bufsp     = ones(Grid.ny, Grid.nx);  % continuity rows
Bound.h_surf    = zeros(1, Grid.nx);       % solid height under each station
Bound.ghost_idx = zeros(1, Grid.nx);       % row of the ghost node, 0 = none
Bound.ghost_wts = zeros(4, Grid.nx);       % its four interpolation weights
% Streamwise stencils next to a face reach into it, as DeHNSSo-v1 did, unless
% Opt.solid_stencil = 'fluid' asks them to stop short.
Bound.stencil_wall = ~(isfield(Opt,'solid_stencil') && strcmpi(Opt.solid_stencil,'fluid'));

if ~strcmpi(Grid.StepType, 'FLAT')
    ny = Grid.ny;   nx = Grid.nx;
    x  = Grid.xun(:).';
    y  = Grid.y(:);                         % wall-normal nodes: row 1 the freestream, row ny the wall
    tol = 1e-9 * max(1, max(abs(y)));       % "on a node" means equal to rounding

    mode = 'auto';
    if isfield(Opt,'solid') && ~isempty(Opt.solid);  mode = lower(Opt.solid);  end
    if strcmp(mode, 'auto')
        snapped = isfield(Grid,'step') && isfield(Grid.step,'snapped') && Grid.step.snapped;
        if snapped;  mode = 'mask';  else;  mode = 'ghost';  end
    end

    % Height of the solid under each station. Each box is widened by the
    % uncertainty gridgen attached to its face, so a station inside that
    % bracket counts as solid, as in DeHNSSo-v1.
    seg = Grid.step.seg;   seg_tol = Grid.step.seg_tol;
    h_surf = zeros(1, nx);
    for b = 1:size(seg,1)
        in = x >= seg(b,1) - seg_tol(b,1) & x <= seg(b,2) + seg_tol(b,1);
        h_surf(in) = max(h_surf(in), seg(b,3));
    end
    Bound.h_surf = h_surf;

    % Station by station. Row 1 is the freestream and row ny the wall, so the
    % first row with y below the surface is the first node inside the solid.
    for i = 1:nx
        h = h_surf(i);
        if h <= 0; continue; end

        k_on    = find(abs(y - h) <= tol, 1);   % a node sits on the surface
        k_below = find(y < h - tol, 1);         % first node inside the solid
        if ~isempty(k_on)
            k_top = k_on;      use_ghost = false;    % no-slip on that node, no ghost row
        else
            k_top = k_below;   use_ghost = strcmp(mode, 'ghost');
        end

        % Velocity rows off from the top node down to the wall. The ghost node is
        % switched off too: LHS replaces its row with the weights below.
        if ~isempty(k_top)
            Bound.bufs(k_top:end, i) = 0;
        end

        % Ghost row: cubic Lagrange through the two fluid nodes above the ghost
        % node, the ghost node and the one below it, evaluated at the surface.
        if use_ghost && ~isempty(k_top) && k_top >= 3 && k_top < ny
            yp = y(k_top-2 : k_top+1);
            w  = zeros(4,1);
            for m = 1:4
                others = yp([1:m-1, m+1:4]);
                w(m)   = prod(h - others) / prod(yp(m) - others);
            end
            Bound.ghost_idx(i)   = k_top;
            Bound.ghost_wts(:,i) = w;
        end

        % Continuity rows off strictly below the surface; the surface node keeps its.
        Bound.bufsp(y < h - tol, i) = 0;
    end

    % At a concave corner, where the solid rises, the continuity row is dropped
    % at the height it rises from: the wall itself for a step standing on it.
    faces  = find(diff(h_surf) ~= 0) + 1;
    rising = faces(h_surf(faces) > h_surf(faces-1));
    for f = rising
        h_up = h_surf(f-1);
        if h_up <= 0
            k = ny;
        else
            k = find(abs(y - h_up) <= tol, 1);
            if isempty(k); [~, k] = min(abs(y - h_up)); end
        end
        Bound.bufsp(k, f) = 0;
    end

    falling = faces(h_surf(faces) < h_surf(faces-1));
    if ~isempty(falling)
        warning('setup_boundary:fallingFace', ...
            ['%d backward-facing face(s) present (first at x = %.4g). The masks are built, ' ...
             'but the concave-corner treatment is applied at rising faces only and this ' ...
             'configuration is not validated.'], numel(falling), x(falling(1)));
    end

    fprintf(' Solid (%s): %s, %d box(es), %d face(s), h = %.4g..%.4g (%d ghost points)\n', ...
        mode, Grid.StepType, size(seg,1), numel(faces), min(seg(:,3)), max(seg(:,3)), nnz(Bound.ghost_idx));
end

%% Inhomogeneous BCs: interpolate Stab.bcw / Stab.bct onto numerical grid
% pchip with extrapolation value 0: the caller's Stab.bcx need not span the
% domain, and stations outside it get zero forcing rather than the NaN that
% 'cubic' without an extrapolation value returns.
is_compr = isfield(Opt,'compressible') && strcmpi(Opt.compressible,'on');
n_bc     = 3 + double(is_compr);   % wall u,v,w[,T]; freestream u,v,w | u,w,rho,T (bc_rows)
nf = size(Stab.bcw, 3);
Bound.bcw = zeros(n_bc, Grid.nx, nf);
Bound.bct = zeros(n_bc, Grid.nx, nf);
has_bcw = any(Stab.bcw(:));
has_bct = any(Stab.bct(:));
if has_bcw || has_bct
    for j = 1:nf
        for k = 1:n_bc
            if has_bcw && any(Stab.bcw(k,:,j))
                Bound.bcw(k,:,j) = interp1(Stab.bcx, Stab.bcw(k,:,j), Grid.xun, 'pchip', 0);
            end
            if has_bct && any(Stab.bct(k,:,j))
                Bound.bct(k,:,j) = interp1(Stab.bcx, Stab.bct(k,:,j), Grid.xun, 'pchip', 0);
            end
        end
    end
    if has_bcw; fprintf(' Inhomogeneous wall BC active (%d modes).\n', nnz(squeeze(any(any(Stab.bcw,1),2)))); end
    if has_bct; fprintf(' Inhomogeneous freestream BC active (%d modes).\n', nnz(squeeze(any(any(Stab.bct,1),2)))); end
end

%% Resonator reactive wall: resolve the admittance and the mouth once
if isfield(Opt, 'hr') && ~isempty(Opt.hr)
    hr = Opt.hr;
    if ~isfield(hr, 'Y') || isempty(hr.Y); hr.Y = hr.mag .* exp(1i*hr.phi); end
    if ~isfield(hr, 'mask') || isempty(hr.mask)
        F       = Grid.forcing(1);
        hr.mask = abs(Grid.xun(1,:) - F.xc) <= F.w;
    end
    hr.mask  = logical(hr.mask(:).');
    Bound.hr = hr;
    fprintf(' Resonator wall BC active: |Y| = %.4f, arg Y = %+.1f deg, mouth = %d stations.\n', ...
            abs(hr.Y(1)), rad2deg(angle(hr.Y(1))), nnz(hr.mask));
end

end

%% Local helpers

function buf = buffer_ramp(buf, ib, nx, kappa)
%BUFFER_RAMP  Apply a tanh ramp over buf(ib:nx), normalized to 1 at ib.
i_vec        = ib:nx;
span         = max(nx - ib, 1);
buf(i_vec)   = 0.5*(1 + tanh(kappa*(1 - 2*(i_vec - ib)/span)));
buf(i_vec)   = buf(i_vec) ./ buf(ib);
end

