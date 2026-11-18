function [Jh, Ja] = nlt_jacobian(q_complex, ctx)
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
% Analytical Jacobian of NLT_HNS with respect to q (2D case only).
%
%   [Jh, Ja] = nlt_jacobian(q_complex, ctx)
%
%   Returns sparse complex (nphys*bsize)^2 matrices holding the Wirtinger
%   derivatives of the NLT forcing,
%
%       d f_NLT = Jh * dq_phys + Ja * conj(dq_phys)
%
%   Jh collecting the terms where the differentiated slot is a physical mode and
%   Ja those where it is the conjugate of one. NB these are the forcing derivatives
%   alone, dr/dQ and dr/dQ*: the paper's J_h and J_a are ML - Jh and -Ja, formed
%   in jacobian_apply. Only the momentum rows (u, v, w) are filled; the inflow
%   rows (station i=1) are zero. Buffer attenuation from ctx.Bound.bufcf applied.
%
%   For nvar=5 the convective momentum forcing is weighted by the base density
%   and keeps the spanwise iu*beta*w_gh coupling:
%
%       xm_j = rho_bar * ( - u_gh*Dx(u_jk) - v_gh*Dy(u_jk) - i*beta_jk*w_gh*u_jk )
%       ym_j = rho_bar * ( - u_gh*Dx(v_jk) - v_gh*Dy(v_jk) - i*beta_jk*w_gh*v_jk )
%       zm_j = rho_bar * ( - u_gh*Dx(w_jk) - v_gh*Dy(w_jk) - i*beta_jk*w_gh*w_jk )
%
%   For nvar=4, rho_bar = 1 and beta = 0, so these reduce to the plain
%   convective bilinear terms.
%
%   ctx fields used: nphys, nf, ny, nx, bsize, phys_modes, HB, dxi, D1, xix,
%   etax, xiy, etay, Bound.bufcf. Optionally cached on first call: ctx.Dx_op,
%   ctx.Dy_op, the one-mode physical derivative operators.

%% Inputs
% Name      size            unit   explanation
% q_complex (nphys*bsize,1) [-]    state of the physical modes as complex coefficients
% ctx       (1,1) struct    [-]    Newton context from newton_context: operators, forcing registry and mode bookkeeping

%% Outputs
% Name size                   unit   explanation
% Jh   (nphys*bsize)^2 sparse [-]    Wirtinger derivative of the nonlinear forcing with respect to q
% Ja   (nphys*bsize)^2 sparse [-]    Wirtinger derivative of the nonlinear forcing with respect to conj(q)

%% Setup: dimensions and mode indexing
nphys = ctx.nphys;
nf    = ctx.nf;
ny    = ctx.ny;
nx    = ctx.nx;
bsize = ctx.bsize;
phys  = ctx.phys_modes;
half  = round(nf/2);
nv    = ctx.nvar;                       % 4 incompressible | 5 compressible
Ntot  = nphys * bsize;

% Grid -> block-local index base for gmap: gbase(k) maps a grid index k in
% 1..ny*nx (column-major) to its offset within a mode block at sub=0, so gmap
% only adds base_off+sub. Avoids a floor/mod on long vectors in the cubic hot path.
gk    = (1:ny*nx).' - 1;
gbase = floor(gk/ny)*nv*ny + mod(gk, ny) + 1;

% Compressible NLT group selector, for debugging and validation: assemble only
% the listed groups (default all). 'a' momentum convection, 'b' continuity,
% 'c' rho'-coupling, 'd' equation-of-state pressure, 'e' viscous variable-mu,
% 'f' energy viscous, 'g' cubic.
jgroups = 'abcdefg';
if isfield(ctx,'jac_groups') && ~isempty(ctx.jac_groups); jgroups = lower(ctx.jac_groups); end
wantg = @(c) any(jgroups==c);
wa = (nv == 4) || wantg('a');   % momentum convection (always on for incompressible)

% Base density for the compressible density-weighted convection rho_bar*(u'.grad)u'.
% Incompressible: rho_bar = 1. Stored column-major (y fastest) to match the
% Dx_op/Dy_op ordering and the ndgrid(1:ny,1:nx) value vectors
if nv == 5
    rhob_vec = ctx.StabGrid.rho(:);
else
    rhob_vec = ones(ny*nx, 1);
end

% Base-flow fields for the rho'-coupling group (c):
% -rho'(d_t + Ubar.grad)q' - rho'(u'.grad)Qbar. The base flow is a quasi-2D
% flat plate, so d_z(base) = 0.
if nv == 5
    G = ctx.StabGrid;
    Ub = G.U(:);  Vb = G.V(:);  Wb = G.W(:);
    bdxU=G.dxU(:); bdyU=G.dyU(:); bdxV=G.dxV(:); bdyV=G.dyV(:);
    bdxW=G.dxW(:); bdyW=G.dyW(:); bdxT=G.dxT(:); bdyT=G.dyT(:);
    z0 = zeros(ny*nx,1);
    bdzU=z0; bdzV=z0; bdzW=z0; bdzT=z0;
    if isfield(G,'dzU')&&~isempty(G.dzU), bdzU=G.dzU(:); end
    if isfield(G,'dzV')&&~isempty(G.dzV), bdzV=G.dzV(:); end
    if isfield(G,'dzW')&&~isempty(G.dzW), bdzW=G.dzW(:); end
    if isfield(G,'dzT')&&~isempty(G.dzT), bdzT=G.dzT(:); end
    % Variable-mu viscous coefficients (groups e/f): mu'(T) = muT*T',
    % lambda = -2/3 mu (Stokes' hypothesis)
    Re_ = G.Re;  RePr_ = G.Re*G.Pr;
    muT_b = G.dmuT(:);  muTT_b = G.d2muT2(:);  mu_b = G.mu(:);
    lamT_b = -2/3*muT_b;  lamTT_b = -2/3*muTT_b;  lam_b = -2/3*mu_b;
    % Physical 2nd-derivative operators, built once
    [Dxx_op, Dyy_op, Dxy_op] = build_diff2_ops(ctx);
    % Spectral index maps for the cubic (3-wave) triples in group (g);
    % mirrors the same indexing in NLT_HNS
    ovc=real(ctx.omegavec(:)).'; bvc=real(ctx.betavec(:)).';
    om0c=min(abs(ovc(abs(ovc)>1e-30))); if isempty(om0c), om0c=1; end
    bt0c=min(abs(bvc(abs(bvc)>1e-30))); if isempty(bt0c), bt0c=1; end
    mmv=round(ovc/om0c); nnv=round(bvc/bt0c);
    mlo=min(mmv); mhi=max(mmv); nlo=min(nnv); nhi=max(nnv);
    Lkp=zeros(mhi-mlo+1, nhi-nlo+1);
    for kk2=1:nf, Lkp(mmv(kk2)-mlo+1, nnv(kk2)-nlo+1)=kk2; end
    RunJg = ctx.Mode.RunJ;
    gamG = ctx.StabGrid.gamma;  MaG = ctx.StabGrid.Ma;
    ikM2g = 1/(gamG*MaG^2);  EcG = (gamG-1)*MaG^2;
end

% One-mode physical derivative operators (size ny*nx, column-major y first)
if isfield(ctx, 'Dx_op') && ~isempty(ctx.Dx_op)
    Dx_op = ctx.Dx_op;
    Dy_op = ctx.Dy_op;
else
    [Dx_op, Dy_op] = build_diff_ops(ctx);
end

% Unpack all physical and conjugate modes into ny x nx arrays (u, v, w[, rho]).
% w is kept even in 2D, since NLT_HNS still computes zm_j, whose derivatives are
% nonzero whenever w ~= 0. rho (compressible) feeds the continuity NLT in group (b).
u_mode = zeros(nf, ny, nx);
v_mode = zeros(nf, ny, nx);
w_mode = zeros(nf, ny, nx);
if nv == 5
    rho_mode = zeros(nf, ny, nx);
    T_mode   = zeros(nf, ny, nx);
end
for k = 1:nphys
    j  = phys(k);
    jc = nf + 1 - j;
    off = (k-1)*bsize;
    qmat = reshape(q_complex(off+1:off+bsize), nv*ny, nx);
    u_mode(j,:,:) = qmat(1:ny,        :);
    v_mode(j,:,:) = qmat(ny+1:2*ny,   :);
    w_mode(j,:,:) = qmat(2*ny+1:3*ny, :);
    if nv == 5
        rho_mode(j,:,:) = qmat(3*ny+1:4*ny, :);
        T_mode(j,:,:)   = qmat(4*ny+1:5*ny, :);
    end
    if jc ~= j
        u_mode(jc,:,:) = conj(qmat(1:ny,        :));
        v_mode(jc,:,:) = conj(qmat(ny+1:2*ny,   :));
        w_mode(jc,:,:) = conj(qmat(2*ny+1:3*ny, :));
        if nv == 5
            rho_mode(jc,:,:) = conj(qmat(3*ny+1:4*ny, :));
            T_mode(jc,:,:)   = conj(qmat(4*ny+1:5*ny, :));
        end
    end
end

% Triplet accumulators, periodically folded into Jh_acc/Ja_acc by fold_tr to
% bound peak memory: the raw list carries ~4x duplicate (row,col) entries, so
% folding sums them and frees the arrays. Triplets are collected in CELL arrays
% and concatenated once per mode, which is O(n); numeric [Ih;rows] growth is
% O(n^2) and used to dominate build time on large M/nx grids.
Ih = {}; Jh_cols = {}; Vh = {};
Ia = {}; Ja_cols = {}; Va = {};
Jh_acc = sparse(Ntot, Ntot);
Ja_acc = sparse(Ntot, Ntot);

% Helper: map mode index ell in [1..nf] to column block in [1..nphys].
% Returns 0 if the mode is not in the active set (skip that entry).
is_conj  = @(ell) (ell <  half);
phys_of  = @(ell) (ell >= half) * ell + (ell < half) * (nf + 1 - ell);

% Build lookup: physical mode index -> column block (0 = inactive)
mode2blk = zeros(nf, 1);
for kk_tmp = 1:nphys
    mode2blk(phys(kk_tmp)) = kk_tmp;
end
col_blk = @(ell) mode2blk(phys_of(ell));

% Operator constants
Zx = sparse(ny*nx, ny*nx);

%% Accumulate the NLT Jacobian per physical mode
for kk = 1:nphys
    j = phys(kk);
    [gharray, jkarray] = find(ctx.HB(:,:,j));
    if isempty(gharray); continue; end

    % Row block offset for this mode's f_j
    row_off = (kk-1)*bsize;
    % u-momentum rows: within this mode's block, the row index of u(y,i) is
    % (i-1)*nv*ny + y (station stride nv*ny; u, v, w sit at sub-offsets
    % 0, ny, 2ny). Maps (y,i) -> global row.
    [Y, X] = ndgrid(1:ny, 1:nx);
    u_row_local = (X(:)-1)*nv*ny +        Y(:);
    v_row_local = (X(:)-1)*nv*ny +   ny + Y(:);
    w_row_local = (X(:)-1)*nv*ny + 2*ny + Y(:);
    u_row = row_off + u_row_local;
    v_row = row_off + v_row_local;
    w_row = row_off + w_row_local;

    % Column-block row mapping: same layout
    u_col_local = u_row_local;
    v_col_local = v_row_local;
    w_col_local = w_row_local;

    % Continuity (rho) + energy (T) rows/cols at sub-offsets 3*ny, 4*ny
    % (compressible only).
    if nv == 5
        rho_row_local = (X(:)-1)*nv*ny + 3*ny + Y(:);
        rho_col_local = rho_row_local;
        rho_row = row_off + rho_row_local;
        T_row_local = (X(:)-1)*nv*ny + 4*ny + Y(:);
        T_col_local = T_row_local;
        T_row = row_off + T_row_local;
    end

    for ip = 1:numel(gharray)
        gh = gharray(ip);
        jk = jkarray(ip);
        beta_jk = ctx.betavec(jk);         % spanwise wavenumber of the jk slot
        if nv == 5; beta_gh = ctx.betavec(gh); end

        % Skip if either mode is not in the active set
        cb_gh_test = col_blk(gh);
        cb_jk_test = col_blk(jk);
        if cb_gh_test == 0 && cb_jk_test == 0; continue; end

        u_gh = squeeze(u_mode(gh,:,:));
        v_gh = squeeze(v_mode(gh,:,:));
        w_gh = squeeze(w_mode(gh,:,:));
        u_jk = squeeze(u_mode(jk,:,:));
        v_jk = squeeze(v_mode(jk,:,:));
        w_jk = squeeze(w_mode(jk,:,:));

        % Physical derivatives of jk-mode fields
        Dx_u_jk = reshape(Dx_op * u_jk(:), ny, nx);
        Dy_u_jk = reshape(Dy_op * u_jk(:), ny, nx);
        Dx_v_jk = reshape(Dx_op * v_jk(:), ny, nx);
        Dy_v_jk = reshape(Dy_op * v_jk(:), ny, nx);
        Dx_w_jk = reshape(Dx_op * w_jk(:), ny, nx);
        Dy_w_jk = reshape(Dy_op * w_jk(:), ny, nx);

        % gh slot, diagonal in (y,x), weighted by rho_bar:
        %   ∂xm_j/∂u_gh = -rho_bar*diag(Dx_u_jk),  ∂xm_j/∂v_gh = -rho_bar*diag(Dy_u_jk)
        %   ∂xm_j/∂w_gh = -rho_bar*i*beta_jk*diag(u_jk)   (compressible only)
        % Rows u/v/w_row, columns u/v/w_col_gh in col_blk(gh).
        cb_gh  = col_blk(gh);
        conj_gh = is_conj(gh);
        if cb_gh > 0 && wa   % gh mode active + group (a) selected
            col_off_gh = (cb_gh-1)*bsize;
            ucol_gh = col_off_gh + u_col_local;
            vcol_gh = col_off_gh + v_col_local;
            wcol_gh = col_off_gh + w_col_local;

            [Ih, Jh_cols, Vh, Ia, Ja_cols, Va] = push_entries( ...
                Ih, Jh_cols, Vh, Ia, Ja_cols, Va, ...
                u_row, ucol_gh, -Dx_u_jk(:).*rhob_vec, conj_gh);
            [Ih, Jh_cols, Vh, Ia, Ja_cols, Va] = push_entries( ...
                Ih, Jh_cols, Vh, Ia, Ja_cols, Va, ...
                u_row, vcol_gh, -Dy_u_jk(:).*rhob_vec, conj_gh);
            [Ih, Jh_cols, Vh, Ia, Ja_cols, Va] = push_entries( ...
                Ih, Jh_cols, Vh, Ia, Ja_cols, Va, ...
                v_row, ucol_gh, -Dx_v_jk(:).*rhob_vec, conj_gh);
            [Ih, Jh_cols, Vh, Ia, Ja_cols, Va] = push_entries( ...
                Ih, Jh_cols, Vh, Ia, Ja_cols, Va, ...
                v_row, vcol_gh, -Dy_v_jk(:).*rhob_vec, conj_gh);
            [Ih, Jh_cols, Vh, Ia, Ja_cols, Va] = push_entries( ...
                Ih, Jh_cols, Vh, Ia, Ja_cols, Va, ...
                w_row, ucol_gh, -Dx_w_jk(:).*rhob_vec, conj_gh);
            [Ih, Jh_cols, Vh, Ia, Ja_cols, Va] = push_entries( ...
                Ih, Jh_cols, Vh, Ia, Ja_cols, Va, ...
                w_row, vcol_gh, -Dy_w_jk(:).*rhob_vec, conj_gh);

            % Spanwise iu*beta_jk*w_gh coupling: ∂(-rho_bar*i*beta*w_gh*q_jk)/∂w_gh
            if beta_jk ~= 0
                cf = -1i*beta_jk .* rhob_vec;
                [Ih, Jh_cols, Vh, Ia, Ja_cols, Va] = push_entries( ...
                    Ih, Jh_cols, Vh, Ia, Ja_cols, Va, ...
                    u_row, wcol_gh, cf.*u_jk(:), conj_gh);
                [Ih, Jh_cols, Vh, Ia, Ja_cols, Va] = push_entries( ...
                    Ih, Jh_cols, Vh, Ia, Ja_cols, Va, ...
                    v_row, wcol_gh, cf.*v_jk(:), conj_gh);
                [Ih, Jh_cols, Vh, Ia, Ja_cols, Va] = push_entries( ...
                    Ih, Jh_cols, Vh, Ia, Ja_cols, Va, ...
                    w_row, wcol_gh, cf.*w_jk(:), conj_gh);
            end
        end

        % jk slot, operator blocks. One operator serves all three momentum rows
        % against their own component (xm<-u_jk, ym<-v_jk, zm<-w_jk):
        %   ∂q_m/∂q_jk = -diag(u_gh)*Dx - diag(v_gh)*Dy - i*beta_jk*diag(w_gh)
        cb_jk   = col_blk(jk);
        conj_jk = is_conj(jk);
        if cb_jk > 0 && wa   % jk mode active + group (a) selected
            col_off_jk = (cb_jk-1)*bsize;

            B_op = -spdiags(u_gh(:), 0, ny*nx, ny*nx) * Dx_op ...
                   -spdiags(v_gh(:), 0, ny*nx, ny*nx) * Dy_op;
            if beta_jk ~= 0
                B_op = B_op - 1i*beta_jk*spdiags(w_gh(:), 0, ny*nx, ny*nx);
            end

            [bi, bj, bv] = find(B_op);
            if ~isempty(bi)
                bv = bv .* rhob_vec(bi);         % rho_bar at the row grid point
                rows_xm = row_off + (floor((bi-1)/ny))*nv*ny + mod(bi-1,ny) + 1;
                cols_uj = col_off_jk + (floor((bj-1)/ny))*nv*ny + mod(bj-1,ny) + 1;
                [Ih, Jh_cols, Vh, Ia, Ja_cols, Va] = push_entries( ...
                    Ih, Jh_cols, Vh, Ia, Ja_cols, Va, ...
                    rows_xm, cols_uj, bv, conj_jk);

                rows_ym = row_off + (floor((bi-1)/ny))*nv*ny + ny + mod(bi-1,ny) + 1;
                cols_vj = col_off_jk + (floor((bj-1)/ny))*nv*ny + ny + mod(bj-1,ny) + 1;
                [Ih, Jh_cols, Vh, Ia, Ja_cols, Va] = push_entries( ...
                    Ih, Jh_cols, Vh, Ia, Ja_cols, Va, ...
                    rows_ym, cols_vj, bv, conj_jk);

                rows_zm = row_off + (floor((bi-1)/ny))*nv*ny + 2*ny + mod(bi-1,ny) + 1;
                cols_wj = col_off_jk + (floor((bj-1)/ny))*nv*ny + 2*ny + mod(bj-1,ny) + 1;
                [Ih, Jh_cols, Vh, Ia, Ja_cols, Va] = push_entries( ...
                    Ih, Jh_cols, Vh, Ia, Ja_cols, Va, ...
                    rows_zm, cols_wj, bv, conj_jk);
            end
        end

        % Shared by groups (b)-(e): gh density, jk temperature, and their gradients.
        if nv == 5
            rho_gh    = squeeze(rho_mode(gh,:,:));
            Dx_rho_gh = reshape(Dx_op*rho_gh(:), ny, nx);   % d rho_gh / dx
            Dy_rho_gh = reshape(Dy_op*rho_gh(:), ny, nx);   % d rho_gh / dy
            T_jk      = squeeze(T_mode(jk,:,:));
            Dx_T_jk   = reshape(Dx_op*T_jk(:), ny, nx);
            Dy_T_jk   = reshape(Dy_op*T_jk(:), ny, nx);
            omega_jk  = ctx.omegavec(jk);
        end

        % Group (b): continuity  c_j = -div(rho' u')
        %   = -[d_x(rho_gh u_jk) + d_y(rho_gh v_jk) + d_z(rho_gh w_jk)]
        % rho-row (+3ny), bilinear in rho_gh (gh) and (u,v,w)_jk (jk).
        if nv == 5 && wantg('b')
            % --- gh slot: d c_j / d rho_gh  (column = rho of gh block) ---
            if cb_gh > 0
                col_off_gh = (cb_gh-1)*bsize;
                Grho = -spdiags(u_jk(:), 0, ny*nx, ny*nx) * Dx_op ...
                       -spdiags(v_jk(:), 0, ny*nx, ny*nx) * Dy_op ...
                       -spdiags(Dx_u_jk(:) + Dy_v_jk(:) ...
                                + 1i*(beta_gh+beta_jk).*w_jk(:), 0, ny*nx, ny*nx);
                [gi, gj, gvv] = find(Grho);
                if ~isempty(gi)
                    rows_c = gmap(gi, gbase, row_off,    3*ny);
                    cols_r = gmap(gj, gbase, col_off_gh, 3*ny);
                    [Ih, Jh_cols, Vh, Ia, Ja_cols, Va] = push_entries( ...
                        Ih, Jh_cols, Vh, Ia, Ja_cols, Va, rows_c, cols_r, gvv, conj_gh);
                end
            end

            % --- jk slot: d c_j / d {u,v,w}_jk  (columns = velocity of jk block) ---
            if cb_jk > 0
                col_off_jk = (cb_jk-1)*bsize;
                % u_jk column: -diag(d rho_gh/dx) - diag(rho_gh)*Dx
                Cu = -spdiags(rho_gh(:),   0, ny*nx, ny*nx) * Dx_op ...
                     -spdiags(Dx_rho_gh(:),0, ny*nx, ny*nx);
                [ci, cj, cvv] = find(Cu);
                if ~isempty(ci)
                    rows_cu = gmap(ci, gbase, row_off,    3*ny);
                    cols_cu = gmap(cj, gbase, col_off_jk, 0);
                    [Ih, Jh_cols, Vh, Ia, Ja_cols, Va] = push_entries( ...
                        Ih, Jh_cols, Vh, Ia, Ja_cols, Va, rows_cu, cols_cu, cvv, conj_jk);
                end
                % v_jk column: -diag(d rho_gh/dy) - diag(rho_gh)*Dy
                Cv = -spdiags(rho_gh(:),   0, ny*nx, ny*nx) * Dy_op ...
                     -spdiags(Dy_rho_gh(:),0, ny*nx, ny*nx);
                [ci, cj, cvv] = find(Cv);
                if ~isempty(ci)
                    rows_cv = gmap(ci, gbase, row_off,    3*ny);
                    cols_cv = gmap(cj, gbase, col_off_jk, ny);
                    [Ih, Jh_cols, Vh, Ia, Ja_cols, Va] = push_entries( ...
                        Ih, Jh_cols, Vh, Ia, Ja_cols, Va, rows_cv, cols_cv, cvv, conj_jk);
                end
                % w_jk column: -i(beta_gh+beta_jk)*diag(rho_gh)  (pure diagonal)
                wcol_jk = col_off_jk + w_col_local;
                [Ih, Jh_cols, Vh, Ia, Ja_cols, Va] = push_entries( ...
                    Ih, Jh_cols, Vh, Ia, Ja_cols, Va, ...
                    rho_row, wcol_jk, -1i*(beta_gh+beta_jk).*rho_gh(:), conj_jk);
            end
        end

        % Group (c): rho'-coupling  q_r = -rho'(d_t q' + Ubar.grad q') - rho'(u'.grad Qbar),
        % d_t -> -i*omega_jk, base d_z = 0. Momentum rows (0, ny, 2ny) and energy (4ny).
        % NOT weighted by rho_bar: rho_gh IS the perturbation density.
        if nv == 5 && wantg('c')
            n2 = ny*nx;
            rvec = rho_gh(:);
            iw = 1i*omega_jk;  ibW = 1i*beta_jk.*Wb;

            % --- gh slot: d q_r / d rho_gh = -diag(L(q_jk)), col = rho of gh ---
            if cb_gh > 0
                col_off_gh = (cb_gh-1)*bsize;
                rhocol_gh = col_off_gh + rho_col_local;
                Bx = -iw.*u_jk(:) + Ub.*Dx_u_jk(:) + Vb.*Dy_u_jk(:) + ibW.*u_jk(:) ...
                     + u_jk(:).*bdxU + v_jk(:).*bdyU + w_jk(:).*bdzU;
                By = -iw.*v_jk(:) + Ub.*Dx_v_jk(:) + Vb.*Dy_v_jk(:) + ibW.*v_jk(:) ...
                     + u_jk(:).*bdxV + v_jk(:).*bdyV + w_jk(:).*bdzV;
                Bz = -iw.*w_jk(:) + Ub.*Dx_w_jk(:) + Vb.*Dy_w_jk(:) + ibW.*w_jk(:) ...
                     + u_jk(:).*bdxW + v_jk(:).*bdyW + w_jk(:).*bdzW;
                BT = -iw.*T_jk(:) + Ub.*Dx_T_jk(:) + Vb.*Dy_T_jk(:) + ibW.*T_jk(:) ...
                     + u_jk(:).*bdxT + v_jk(:).*bdyT + w_jk(:).*bdzT;
                [Ih,Jh_cols,Vh,Ia,Ja_cols,Va]=push_entries(Ih,Jh_cols,Vh,Ia,Ja_cols,Va, u_row, rhocol_gh, -Bx, conj_gh);
                [Ih,Jh_cols,Vh,Ia,Ja_cols,Va]=push_entries(Ih,Jh_cols,Vh,Ia,Ja_cols,Va, v_row, rhocol_gh, -By, conj_gh);
                [Ih,Jh_cols,Vh,Ia,Ja_cols,Va]=push_entries(Ih,Jh_cols,Vh,Ia,Ja_cols,Va, w_row, rhocol_gh, -Bz, conj_gh);
                [Ih,Jh_cols,Vh,Ia,Ja_cols,Va]=push_entries(Ih,Jh_cols,Vh,Ia,Ja_cols,Va, T_row, rhocol_gh, -BT, conj_gh);
            end

            % --- jk slot: d q_r / d q_jk = -diag(rho_gh)*L ---
            if cb_jk > 0
                col_off_jk = (cb_jk-1)*bsize;
                ucol_jk = col_off_jk + u_col_local;
                vcol_jk = col_off_jk + v_col_local;
                wcol_jk = col_off_jk + w_col_local;
                Rg = spdiags(rvec, 0, n2, n2);
                % common base-convective operator: -i*omega + Ubar Dx + Vbar Dy + i*beta*Wbar
                Cc = -iw*speye(n2) + spdiags(Ub,0,n2,n2)*Dx_op ...
                     + spdiags(Vb,0,n2,n2)*Dy_op + spdiags(ibW,0,n2,n2);
                % self operators (+ base-gradient self diagonal): x<-u, y<-v, z<-w, E<-T
                selfM = { -Rg*(Cc + spdiags(bdxU,0,n2,n2)), 0; ...
                          -Rg*(Cc + spdiags(bdyV,0,n2,n2)), ny; ...
                          -Rg*(Cc + spdiags(bdzW,0,n2,n2)), 2*ny; ...
                          -Rg*Cc,                            4*ny };
                for s = 1:4
                    [oi,oj,ov] = find(selfM{s,1});
                    if isempty(oi); continue; end
                    sub = selfM{s,2};
                    [Ih,Jh_cols,Vh,Ia,Ja_cols,Va]=push_entries(Ih,Jh_cols,Vh,Ia,Ja_cols,Va, ...
                        gmap(oi,gbase,row_off,sub), gmap(oj,gbase,col_off_jk,sub), ov, conj_jk);
                end
                % cross base-gradient diagonals (row <- other-field col)
                [Ih,Jh_cols,Vh,Ia,Ja_cols,Va]=push_entries(Ih,Jh_cols,Vh,Ia,Ja_cols,Va, u_row, vcol_jk, -rvec.*bdyU, conj_jk);
                [Ih,Jh_cols,Vh,Ia,Ja_cols,Va]=push_entries(Ih,Jh_cols,Vh,Ia,Ja_cols,Va, u_row, wcol_jk, -rvec.*bdzU, conj_jk);
                [Ih,Jh_cols,Vh,Ia,Ja_cols,Va]=push_entries(Ih,Jh_cols,Vh,Ia,Ja_cols,Va, v_row, ucol_jk, -rvec.*bdxV, conj_jk);
                [Ih,Jh_cols,Vh,Ia,Ja_cols,Va]=push_entries(Ih,Jh_cols,Vh,Ia,Ja_cols,Va, v_row, wcol_jk, -rvec.*bdzV, conj_jk);
                [Ih,Jh_cols,Vh,Ia,Ja_cols,Va]=push_entries(Ih,Jh_cols,Vh,Ia,Ja_cols,Va, w_row, ucol_jk, -rvec.*bdxW, conj_jk);
                [Ih,Jh_cols,Vh,Ia,Ja_cols,Va]=push_entries(Ih,Jh_cols,Vh,Ia,Ja_cols,Va, w_row, vcol_jk, -rvec.*bdyW, conj_jk);
                [Ih,Jh_cols,Vh,Ia,Ja_cols,Va]=push_entries(Ih,Jh_cols,Vh,Ia,Ja_cols,Va, T_row, ucol_jk, -rvec.*bdxT, conj_jk);
                [Ih,Jh_cols,Vh,Ia,Ja_cols,Va]=push_entries(Ih,Jh_cols,Vh,Ia,Ja_cols,Va, T_row, vcol_jk, -rvec.*bdyT, conj_jk);
                [Ih,Jh_cols,Vh,Ia,Ja_cols,Va]=push_entries(Ih,Jh_cols,Vh,Ia,Ja_cols,Va, T_row, wcol_jk, -rvec.*bdzT, conj_jk);
            end
        end

        % Group (d): nonlinear EOS pressure. p' is linear in (rho',T'), so the
        % quadratic rho'T' piece forces: momentum -grad(rho'T')/(gamma Ma^2),
        % energy +Ec*(d_t+Ubar.grad)(rho'T')/(gamma Ma^2), d_t -> -i(om_gh+om_jk).
        % Bilinear in rho_gh (gh, rho-col, from group b) and T_jk (jk, T-col, from c).
        if nv == 5 && wantg('d')
            n2   = ny*nx;
            gam  = ctx.StabGrid.gamma;  Ma = ctx.StabGrid.Ma;
            ikM2 = 1/(gam*Ma^2);        Ec = (gam-1)*Ma^2;
            om_gh = ctx.omegavec(gh);   om_jk = ctx.omegavec(jk);
            bsum  = beta_gh + beta_jk;
            Tv  = T_jk(:);      rgv = rho_gh(:);
            DxT = Dx_T_jk(:);   DyT = Dy_T_jk(:);
            Dxr = Dx_rho_gh(:);  Dyr = Dy_rho_gh(:);

            % --- gh slot: d/d rho_gh  (column = rho of gh block) ---
            if cb_gh > 0
                col_off_gh = (cb_gh-1)*bsize;
                rhocol_gh = col_off_gh + rho_col_local;
                Ox = -ikM2*( spdiags(Tv,0,n2,n2)*Dx_op + spdiags(DxT,0,n2,n2) );
                Oy = -ikM2*( spdiags(Tv,0,n2,n2)*Dy_op + spdiags(DyT,0,n2,n2) );
                Eg = Ec*ikM2*( spdiags(-1i*(om_gh+om_jk).*Tv,0,n2,n2) ...
                     + spdiags(Ub.*Tv,0,n2,n2)*Dx_op + spdiags(Ub.*DxT,0,n2,n2) ...
                     + spdiags(Vb.*Tv,0,n2,n2)*Dy_op + spdiags(Vb.*DyT,0,n2,n2) ...
                     + spdiags(Wb.*(1i*bsum).*Tv,0,n2,n2) );
                gdM = { Ox,0; Oy,ny; Eg,4*ny };
                for s = 1:size(gdM,1)
                    [oi,oj,ov] = find(gdM{s,1});
                    if isempty(oi); continue; end
                    [Ih,Jh_cols,Vh,Ia,Ja_cols,Va]=push_entries(Ih,Jh_cols,Vh,Ia,Ja_cols,Va, ...
                        gmap(oi,gbase,row_off,gdM{s,2}), gmap(oj,gbase,col_off_gh,3*ny), ov, conj_gh);
                end
                % z-mom diagonal: -ikM2*i(beta_gh+beta_jk)*diag(T_jk)
                [Ih,Jh_cols,Vh,Ia,Ja_cols,Va]=push_entries(Ih,Jh_cols,Vh,Ia,Ja_cols,Va, ...
                    w_row, rhocol_gh, -ikM2*1i*bsum.*Tv, conj_gh);
            end

            % --- jk slot: d/d T_jk  (column = T of jk block) ---
            if cb_jk > 0
                col_off_jk = (cb_jk-1)*bsize;
                Tcol_jk = col_off_jk + T_col_local;
                Px = -ikM2*( spdiags(Dxr,0,n2,n2) + spdiags(rgv,0,n2,n2)*Dx_op );
                Py = -ikM2*( spdiags(Dyr,0,n2,n2) + spdiags(rgv,0,n2,n2)*Dy_op );
                Ej = Ec*ikM2*( spdiags(-1i*(om_gh+om_jk).*rgv,0,n2,n2) ...
                     + spdiags(Ub.*Dxr,0,n2,n2) + spdiags(Ub.*rgv,0,n2,n2)*Dx_op ...
                     + spdiags(Vb.*Dyr,0,n2,n2) + spdiags(Vb.*rgv,0,n2,n2)*Dy_op ...
                     + spdiags(Wb.*(1i*bsum).*rgv,0,n2,n2) );
                gdN = { Px,0; Py,ny; Ej,4*ny };
                for s = 1:size(gdN,1)
                    [oi,oj,ov] = find(gdN{s,1});
                    if isempty(oi); continue; end
                    [Ih,Jh_cols,Vh,Ia,Ja_cols,Va]=push_entries(Ih,Jh_cols,Vh,Ia,Ja_cols,Va, ...
                        gmap(oi,gbase,row_off,gdN{s,2}), gmap(oj,gbase,col_off_jk,4*ny), ov, conj_jk);
                end
                % z-mom diagonal: -ikM2*i(beta_gh+beta_jk)*diag(rho_gh)
                [Ih,Jh_cols,Vh,Ia,Ja_cols,Va]=push_entries(Ih,Jh_cols,Vh,Ia,Ja_cols,Va, ...
                    w_row, Tcol_jk, -ikM2*1i*bsum.*rgv, conj_jk);
            end
        end

        % Group (e): variable-mu viscous momentum forcing (NLT_HNS 226-280). Bilinear
        % in mu'=muT*T_gh
        % (gh, T-col) and the velocity_jk 1st/2nd derivatives (jk, u/v/w-cols). Each
        % term is C(T_gh).*(op*vel_jk): d/dT_gh gives diag/Dx/Dy/i*beta_gh, d/dvel
        % gives diag(C)*op.
        if nv == 5 && wantg('e')
            n2 = ny*nx;
            dg = @(f) spdiags(f(:),0,n2,n2);
            SDx=Dx_op; SDy=Dy_op; SDxx=Dxx_op; SDyy=Dyy_op; SDxy=Dxy_op; In=speye(n2);
            b = beta_jk; bg = beta_gh; b2 = b^2; reI = 1/Re_;
            T_ghf = squeeze(T_mode(gh,:,:));  Tg = T_ghf(:);
            dTgx = Dx_op*Tg;  dTgy = Dy_op*Tg;
            uj=u_jk(:); vj=v_jk(:); wj=w_jk(:);
            ux=Dx_u_jk(:); uy=Dy_u_jk(:); vx=Dx_v_jk(:); vy=Dy_v_jk(:); wx=Dx_w_jk(:); wy=Dy_w_jk(:);
            uxx=Dxx_op*uj; uyy=Dyy_op*uj; uxy=Dxy_op*uj;
            vxx=Dxx_op*vj; vyy=Dyy_op*vj; vxy=Dxy_op*vj;
            wxx=Dxx_op*wj; wyy=Dyy_op*wj;
            divU = ux + vy + 1i*b.*wj;

            % ---- x-momentum ----
            cT1=reI*lamT_b.*Tg; cT2=2*reI*muT_b.*Tg; cT3=reI*muT_b.*Tg;
            cT4=reI*(lamTT_b.*bdxT.*Tg + lamT_b.*dTgx);
            cT5=2*reI*(muTT_b.*bdxT.*Tg + muT_b.*dTgx);
            cT6=reI*(muTT_b.*bdyT.*Tg + muT_b.*dTgy);
            cT7=reI*(muTT_b.*bdzT.*Tg + muT_b.*(1i*bg).*Tg);
            f1=uxx+vxy+1i*b.*wx; f2=uxx; f3=uyy+vxy-b2*uj+1i*b.*wx;
            f4=divU; f5=ux; f6=uy+vx; f7=1i*b.*uj+wx;
            J_xu = dg(cT1+cT2)*SDxx + dg(cT3)*(SDyy-b2*In) + dg(cT4+cT5)*SDx + dg(cT6)*SDy + (1i*b)*dg(cT7);
            J_xv = dg(cT1+cT3)*SDxy + dg(cT4)*SDy + dg(cT6)*SDx;
            J_xw = (1i*b)*dg(cT1+cT3)*SDx + (1i*b)*dg(cT4) + dg(cT7)*SDx;
            Gx_d  = reI*lamT_b.*f1 + 2*reI*muT_b.*f2 + reI*muT_b.*f3 ...
                  + reI*lamTT_b.*bdxT.*f4 + 2*reI*muTT_b.*bdxT.*f5 + reI*muTT_b.*bdyT.*f6 + reI*muTT_b.*bdzT.*f7;
            Gx_dx = reI*lamT_b.*f4 + 2*reI*muT_b.*f5;
            Gx_dy = reI*muT_b.*f6;
            Gx_bg = reI*muT_b.*f7;
            G_x = dg(Gx_d) + dg(Gx_dx)*SDx + dg(Gx_dy)*SDy + (1i*bg)*dg(Gx_bg);

            % ---- y-momentum ----
            cS1=reI*lamT_b.*Tg; cS2=reI*muT_b.*Tg; cS3=reI*muT_b.*Tg;
            cS4=reI*(lamTT_b.*bdyT.*Tg + lamT_b.*dTgy);
            cS5=2*reI*(muTT_b.*bdyT.*Tg + muT_b.*dTgy);
            cS6=reI*(muTT_b.*bdxT.*Tg + muT_b.*dTgx);
            cS7=reI*(muTT_b.*bdzT.*Tg + muT_b.*(1i*bg).*Tg);
            g1=uxy+vyy+1i*b.*wy; g2=2*vyy+uxy; g3=vxx-b2*vj+1i*b.*wy;
            g4=divU; g5=vy; g6=uy+vx; g7=1i*b.*vj+wy;
            J_yu = dg(cS1+cS2)*SDxy + dg(cS6)*SDy + dg(cS4)*SDx;  % cS4*SDx from divU d/du
            J_yv = dg(cS1+2*cS2)*SDyy + dg(cS3)*(SDxx-b2*In) + dg(cS4+cS5)*SDy + dg(cS6)*SDx + (1i*b)*dg(cS7);
            J_yw = (1i*b)*dg(cS1+cS3)*SDy + (1i*b)*dg(cS4) + dg(cS7)*SDy;
            Gy_d  = reI*lamT_b.*g1 + reI*muT_b.*g2 + reI*muT_b.*g3 ...
                  + reI*lamTT_b.*bdyT.*g4 + 2*reI*muTT_b.*bdyT.*g5 + reI*muTT_b.*bdxT.*g6 + reI*muTT_b.*bdzT.*g7;
            Gy_dx = reI*muT_b.*g6;
            Gy_dy = reI*lamT_b.*g4 + 2*reI*muT_b.*g5;
            Gy_bg = reI*muT_b.*g7;
            G_y = dg(Gy_d) + dg(Gy_dx)*SDx + dg(Gy_dy)*SDy + (1i*bg)*dg(Gy_bg);

            % ---- z-momentum ----
            cR1=reI*lamT_b.*Tg; cR2=reI*muT_b.*Tg; cR3=reI*muT_b.*Tg;
            cR4=reI*(lamTT_b.*bdzT.*Tg + lamT_b.*(1i*bg).*Tg);
            cR5=2*reI*(muTT_b.*bdzT.*Tg + muT_b.*(1i*bg).*Tg);
            cR6=reI*(muTT_b.*bdxT.*Tg + muT_b.*dTgx);
            cR7=reI*(muTT_b.*bdyT.*Tg + muT_b.*dTgy);
            h1=1i*b.*ux+1i*b.*vy-b2*wj; h2=-2*b2*wj+1i*b.*ux; h3=wxx+1i*b.*vy+wyy;
            h4=divU; h5=1i*b.*wj; h6=1i*b.*uj+wx; h7=1i*b.*vj+wy;
            J_zu = (1i*b)*dg(cR1+cR2)*SDx + dg(cR4)*SDx + (1i*b)*dg(cR6);
            J_zv = (1i*b)*dg(cR1+cR3)*SDy + dg(cR4)*SDy + (1i*b)*dg(cR7);
            J_zw = dg(cR1)*(-b2*In) + dg(cR2)*(-2*b2*In) + dg(cR3)*(SDxx+SDyy) ...
                 + (1i*b)*dg(cR4+cR5) + dg(cR6)*SDx + dg(cR7)*SDy;
            Gz_d  = reI*lamT_b.*h1 + reI*muT_b.*h2 + reI*muT_b.*h3 ...
                  + reI*lamTT_b.*bdzT.*h4 + 2*reI*muTT_b.*bdzT.*h5 + reI*muTT_b.*bdxT.*h6 + reI*muTT_b.*bdyT.*h7;
            Gz_dx = reI*muT_b.*h6;
            Gz_dy = reI*muT_b.*h7;
            Gz_bg = reI*lamT_b.*h4 + 2*reI*muT_b.*h5;
            G_z = dg(Gz_d) + dg(Gz_dx)*SDx + dg(Gz_dy)*SDy + (1i*bg)*dg(Gz_bg);

            % ---- push gh-slot (T-col of gh block) ----
            if cb_gh > 0
                col_off_gh = (cb_gh-1)*bsize;
                geG = { G_x,0; G_y,ny; G_z,2*ny };
                for s = 1:3
                    [oi,oj,ov] = find(geG{s,1});
                    if isempty(oi); continue; end
                    [Ih,Jh_cols,Vh,Ia,Ja_cols,Va]=push_entries(Ih,Jh_cols,Vh,Ia,Ja_cols,Va, ...
                        gmap(oi,gbase,row_off,geG{s,2}), gmap(oj,gbase,col_off_gh,4*ny), ov, conj_gh);
                end
            end
            % ---- push jk-slot (u/v/w-col of jk block) ----
            if cb_jk > 0
                col_off_jk = (cb_jk-1)*bsize;
                geJ = { J_xu,0,0; J_xv,0,ny; J_xw,0,2*ny; ...
                        J_yu,ny,0; J_yv,ny,ny; J_yw,ny,2*ny; ...
                        J_zu,2*ny,0; J_zv,2*ny,ny; J_zw,2*ny,2*ny };
                for s = 1:9
                    [oi,oj,ov] = find(geJ{s,1});
                    if isempty(oi); continue; end
                    [Ih,Jh_cols,Vh,Ia,Ja_cols,Va]=push_entries(Ih,Jh_cols,Vh,Ia,Ja_cols,Va, ...
                        gmap(oi,gbase,row_off,geJ{s,2}), gmap(oj,gbase,col_off_jk,geJ{s,3}), ov, conj_jk);
                end
            end
        end

        % Group (f): energy viscous forcing (NLT_HNS 282-302), T-row only. Heat
        % conduction (k=mu), bilinear in k'=muT*T_gh (gh, T-col) and the T_jk 2nd
        % derivatives (jk, T-col), factor 1/(Re*Pr).
        if nv == 5 && wantg('f')
            n2 = ny*nx;
            dg = @(f) spdiags(f(:),0,n2,n2);
            SDx=Dx_op; SDy=Dy_op; SDxx=Dxx_op; SDyy=Dyy_op; In=speye(n2);
            b=beta_jk; bg=beta_gh; b2=b^2; rp=1/RePr_;
            T_ghf = squeeze(T_mode(gh,:,:)); Tg=T_ghf(:);
            dTgx = Dx_op*Tg;  dTgy = Dy_op*Tg;
            Tv=T_jk(:); DxT=Dx_T_jk(:); DyT=Dy_T_jk(:);
            Txx=Dxx_op*Tv; Tyy=Dyy_op*Tv;
            % --- conduction ---
            cC1 = rp*(muTT_b.*bdxT.*Tg + muT_b.*dTgx);   % coeff of dT_jk_dx
            cC3 = rp*(muTT_b.*bdyT.*Tg + muT_b.*dTgy);   % coeff of dT_jk_dy
            cCm = rp*muT_b.*Tg;                           % coeff of Txx, Tyy, -b2*T
            cC5 = rp*(muTT_b.*bdzT.*Tg + muT_b.*(1i*bg).*Tg);   % coeff of i*b*T_jk
            % jk-slot d/dT_jk (T-col)
            J_ET = dg(cC1)*SDx + dg(cC3)*SDy + dg(cCm)*(SDxx + SDyy - b2*In) + (1i*b)*dg(cC5);
            % gh-slot d/dT_gh (T-col): fields multiplying T_gh-linear coefficients
            fC1=DxT; fC2=Txx; fC3=DyT; fC4=Tyy; fC5=1i*b.*Tv; fC6=-b2.*Tv;
            GET_d = rp*muTT_b.*bdxT.*fC1 + rp*muT_b.*fC2 + rp*muTT_b.*bdyT.*fC3 ...
                  + rp*muT_b.*fC4 + rp*muTT_b.*bdzT.*fC5 + rp*muT_b.*fC6;
            GET_dx = rp*muT_b.*fC1;   GET_dy = rp*muT_b.*fC3;   GET_bg = rp*muT_b.*fC5;
            G_ET = dg(GET_d) + dg(GET_dx)*SDx + dg(GET_dy)*SDy + (1i*bg)*dg(GET_bg);

            % --- dissipation D-B: muT*T_gh * base-grad * jk-vel-grad ---
            gam=ctx.StabGrid.gamma; Ma=ctx.StabGrid.Ma; Ec=(gam-1)*Ma^2; er=Ec/Re_;
            ux=Dx_u_jk(:); uy=Dy_u_jk(:); vx=Dx_v_jk(:); vy=Dy_v_jk(:); wx=Dx_w_jk(:); wy=Dy_w_jk(:);
            uj=u_jk(:); vj=v_jk(:); wj=w_jk(:);
            divU = ux + vy + 1i*b.*wj;
            dilB = bdxU + bdyV + bdzW;                    % base dilatation
            GDB = er*lamT_b.*(2*dilB).*divU ...
                + er*4*muT_b.*(bdxU.*ux + bdyV.*vy + bdzW.*(1i*b.*wj)) ...
                + er*2*muT_b.*((bdyU+bdxV).*(uy+vx)) ...
                + er*2*muT_b.*(bdyW.*wy + bdyW.*(1i*b.*vj)) ...
                + er*2*muT_b.*(bdxW.*(1i*b.*wj + wx) + bdxW.*(1i*b.*uj));
            G_ET = G_ET + dg(GDB);
            cb1 = er*lamT_b.*(2*dilB).*Tg;                % B1 divU coeff (T_gh-weighted)
            cbM = er*muT_b.*Tg;                           % base for B2..B5
            J_Eu = dg(cb1)*SDx + dg(4*cbM.*bdxU)*SDx + dg(2*cbM.*(bdyU+bdxV))*SDy + (1i*b)*dg(2*cbM.*bdxW);
            J_Ev = dg(cb1)*SDy + dg(4*cbM.*bdyV)*SDy + dg(2*cbM.*(bdyU+bdxV))*SDx + (1i*b)*dg(2*cbM.*bdyW);
            J_Ew = (1i*b)*dg(cb1) + (1i*b)*dg(4*cbM.*bdzW) + dg(2*cbM.*bdyW)*SDy ...
                 + dg(2*cbM.*bdxW)*(SDx + 1i*b*In);

            % --- dissipation D-A: (lam,mu)_b * grad(vel_gh) * grad(vel_jk) ---
            ugv=u_gh(:); vgv=v_gh(:); wgv=w_gh(:);
            dugx=Dx_op*ugv; dugy=Dy_op*ugv; dvgx=Dx_op*vgv; dvgy=Dy_op*vgv; dwgx=Dx_op*wgv; dwgy=Dy_op*wgv;
            ibg=1i*bg; ib=1i*b;
            % gh-slot vel-col operators (E-row <- u/v/w_gh)
            G_Eug = dg(er*lam_b.*(ux+2*vy+2*ib.*wj) + 2*er*mu_b.*ux)*SDx + dg(er*mu_b.*(uy+2*vx))*SDy;
            G_Evg = dg(er*lam_b.*(vy+2*ib.*wj) + 2*er*mu_b.*vy)*SDy + dg(er*mu_b.*vx)*SDx ...
                  + ibg*dg(er*mu_b.*(ib.*vj));
            G_Ewg = ibg*dg(er*lam_b.*(ib.*wj) + 2*er*mu_b.*(ib.*wj)) ...
                  + dg(er*mu_b.*(wy+2*ib.*vj))*SDy + dg(er*mu_b.*(wx+2*ib.*uj))*SDx;
            % jk-slot vel-col operators (add D-A to D-B J_Eu/J_Ev/J_Ew)
            J_Eu = J_Eu + dg(er*(lam_b+2*mu_b).*dugx)*SDx + dg(er*mu_b.*dugy)*SDy ...
                 + ib*dg(2*er*mu_b.*dwgx + er*mu_b.*(ibg.*ugv));
            J_Ev = J_Ev + dg(er*lam_b.*(dvgy+2*dugx) + 2*er*mu_b.*dvgy)*SDy ...
                 + dg(er*mu_b.*(2*dugy+dvgx))*SDx + ib*dg(2*er*mu_b.*dwgy + er*mu_b.*(ibg.*vgv));
            J_Ew = J_Ew + dg(er*mu_b.*dwgy)*SDy + dg(er*mu_b.*dwgx)*SDx ...
                 + ib*dg(er*lam_b.*(ibg.*wgv) + 2*er*lam_b.*(dugx+dvgy) + 2*er*mu_b.*(ibg.*wgv));

            % push (energy row = 4ny)
            if cb_gh > 0
                col_off_gh = (cb_gh-1)*bsize;
                gfG = { G_ET,4*ny; G_Eug,0; G_Evg,ny; G_Ewg,2*ny };
                for s = 1:4
                    [oi,oj,ov]=find(gfG{s,1});
                    if isempty(oi); continue; end
                    [Ih,Jh_cols,Vh,Ia,Ja_cols,Va]=push_entries(Ih,Jh_cols,Vh,Ia,Ja_cols,Va, ...
                        gmap(oi,gbase,row_off,4*ny), gmap(oj,gbase,col_off_gh,gfG{s,2}), ov, conj_gh);
                end
            end
            if cb_jk > 0
                col_off_jk = (cb_jk-1)*bsize;
                gfJ = { J_ET,4*ny; J_Eu,0; J_Ev,ny; J_Ew,2*ny };
                for s = 1:4
                    [oi,oj,ov]=find(gfJ{s,1});
                    if isempty(oi); continue; end
                    [Ih,Jh_cols,Vh,Ia,Ja_cols,Va]=push_entries(Ih,Jh_cols,Vh,Ia,Ja_cols,Va, ...
                        gmap(oi,gbase,row_off,4*ny), gmap(oj,gbase,col_off_jk,gfJ{s,2}), ov, conj_jk);
                end
            end
        end
    end

    % Group (g): cubic term, trilinear, summed over triples (a,b,c) whose spectral
    % indices add to j. First sub-part -rho_a (u_b.grad) q_c writes momentum and
    % energy rows; slots a (rho-col), b (velocity-cols), c (q-col, operator).
    if nv == 5 && wantg('g')
        dg = @(f) spdiags(f(:),0,ny*nx,ny*nx);
        j = phys(kk);  mj=mmv(j); nj=nnv(j);
        for a = RunJg
          cba=col_blk(a); conj_a=is_conj(a);
          for b = RunJg
            m3i=mj-mmv(a)-mmv(b); n3i=nj-nnv(a)-nnv(b);
            if m3i<mlo||m3i>mhi||n3i<nlo||n3i>nhi, continue; end
            c=Lkp(m3i-mlo+1, n3i-nlo+1);
            if c==0 || ~any(RunJg==c), continue; end
            cbb=col_blk(b); conj_b=is_conj(b);
            cbc=col_blk(c); conj_c=is_conj(c);
            r1 = squeeze(rho_mode(a,:,:)); r1=r1(:);
            u2=squeeze(u_mode(b,:,:)); v2=squeeze(v_mode(b,:,:)); w2=squeeze(w_mode(b,:,:));
            u2v=u2(:); v2v=v2(:); w2v=w2(:);
            u3=squeeze(u_mode(c,:,:)); v3=squeeze(v_mode(c,:,:)); w3=squeeze(w_mode(c,:,:)); T3=squeeze(T_mode(c,:,:));
            u3v=u3(:); v3v=v3(:); w3v=w3(:); T3v=T3(:);
            bc=ctx.betavec(c); ibc=1i*bc;
            du3x=Dx_op*u3v; du3y=Dy_op*u3v; dv3x=Dx_op*v3v; dv3y=Dy_op*v3v;
            dw3x=Dx_op*w3v; dw3y=Dy_op*w3v; dT3x=Dx_op*T3v; dT3y=Dy_op*T3v;

            % c-slot: shared convective operator -rho1*(u2*Dx + v2*Dy + ibc*w2)
            if cbc > 0
                col_off_c=(cbc-1)*bsize;
                Op_c = -dg(r1.*u2v)*Dx_op - dg(r1.*v2v)*Dy_op - dg(r1.*(ibc.*w2v));
                [oi,oj,ov]=find(Op_c);
                if ~isempty(oi)
                    cc = [0 0; ny ny; 2*ny 2*ny; 4*ny 4*ny];  % [row_sub col_sub]: x<-u,y<-v,z<-w,E<-T
                    for s=1:4
                        [Ih,Jh_cols,Vh,Ia,Ja_cols,Va]=push_entries(Ih,Jh_cols,Vh,Ia,Ja_cols,Va, ...
                            gmap(oi,gbase,row_off,cc(s,1)), gmap(oj,gbase,col_off_c,cc(s,2)), ov, conj_c);
                    end
                end
            end
            % a-slot: d/d rho_a = -(bracket), diagonal into rho-col of a
            if cba > 0
                rhocol_a=(cba-1)*bsize + rho_col_local;
                brx=-(u2v.*du3x+v2v.*du3y+ibc.*w2v.*u3v);
                bry=-(u2v.*dv3x+v2v.*dv3y+ibc.*w2v.*v3v);
                brz=-(u2v.*dw3x+v2v.*dw3y+ibc.*w2v.*w3v);
                bre=-(u2v.*dT3x+v2v.*dT3y+ibc.*w2v.*T3v);
                [Ih,Jh_cols,Vh,Ia,Ja_cols,Va]=push_entries(Ih,Jh_cols,Vh,Ia,Ja_cols,Va, u_row, rhocol_a, brx, conj_a);
                [Ih,Jh_cols,Vh,Ia,Ja_cols,Va]=push_entries(Ih,Jh_cols,Vh,Ia,Ja_cols,Va, v_row, rhocol_a, bry, conj_a);
                [Ih,Jh_cols,Vh,Ia,Ja_cols,Va]=push_entries(Ih,Jh_cols,Vh,Ia,Ja_cols,Va, w_row, rhocol_a, brz, conj_a);
                [Ih,Jh_cols,Vh,Ia,Ja_cols,Va]=push_entries(Ih,Jh_cols,Vh,Ia,Ja_cols,Va, T_row, rhocol_a, bre, conj_a);
            end
            % b-slot: d/d u_b,v_b,w_b = -rho1*(grad q3 components), diagonals
            if cbb > 0
                col_off_b=(cbb-1)*bsize;
                ucol_b=col_off_b+u_col_local; vcol_b=col_off_b+v_col_local; wcol_b=col_off_b+w_col_local;
                % x-row
                [Ih,Jh_cols,Vh,Ia,Ja_cols,Va]=push_entries(Ih,Jh_cols,Vh,Ia,Ja_cols,Va, u_row, ucol_b, -r1.*du3x, conj_b);
                [Ih,Jh_cols,Vh,Ia,Ja_cols,Va]=push_entries(Ih,Jh_cols,Vh,Ia,Ja_cols,Va, u_row, vcol_b, -r1.*du3y, conj_b);
                [Ih,Jh_cols,Vh,Ia,Ja_cols,Va]=push_entries(Ih,Jh_cols,Vh,Ia,Ja_cols,Va, u_row, wcol_b, -r1.*(ibc.*u3v), conj_b);
                % y-row
                [Ih,Jh_cols,Vh,Ia,Ja_cols,Va]=push_entries(Ih,Jh_cols,Vh,Ia,Ja_cols,Va, v_row, ucol_b, -r1.*dv3x, conj_b);
                [Ih,Jh_cols,Vh,Ia,Ja_cols,Va]=push_entries(Ih,Jh_cols,Vh,Ia,Ja_cols,Va, v_row, vcol_b, -r1.*dv3y, conj_b);
                [Ih,Jh_cols,Vh,Ia,Ja_cols,Va]=push_entries(Ih,Jh_cols,Vh,Ia,Ja_cols,Va, v_row, wcol_b, -r1.*(ibc.*v3v), conj_b);
                % z-row
                [Ih,Jh_cols,Vh,Ia,Ja_cols,Va]=push_entries(Ih,Jh_cols,Vh,Ia,Ja_cols,Va, w_row, ucol_b, -r1.*dw3x, conj_b);
                [Ih,Jh_cols,Vh,Ia,Ja_cols,Va]=push_entries(Ih,Jh_cols,Vh,Ia,Ja_cols,Va, w_row, vcol_b, -r1.*dw3y, conj_b);
                [Ih,Jh_cols,Vh,Ia,Ja_cols,Va]=push_entries(Ih,Jh_cols,Vh,Ia,Ja_cols,Va, w_row, wcol_b, -r1.*(ibc.*w3v), conj_b);
                % E-row
                [Ih,Jh_cols,Vh,Ia,Ja_cols,Va]=push_entries(Ih,Jh_cols,Vh,Ia,Ja_cols,Va, T_row, ucol_b, -r1.*dT3x, conj_b);
                [Ih,Jh_cols,Vh,Ia,Ja_cols,Va]=push_entries(Ih,Jh_cols,Vh,Ia,Ja_cols,Va, T_row, vcol_b, -r1.*dT3y, conj_b);
                [Ih,Jh_cols,Vh,Ia,Ja_cols,Va]=push_entries(Ih,Jh_cols,Vh,Ia,Ja_cols,Va, T_row, wcol_b, -r1.*(ibc.*T3v), conj_b);
            end

            % --- cubic dissipation: (lamT/muT)*T_a * strain_b * strain_c (E-row) ---
            T1cv=reshape(T_mode(a,:,:),ny,nx); T1cv=T1cv(:);
            bb=ctx.betavec(b);
            du2x=Dx_op*u2v; du2y=Dy_op*u2v; dv2x=Dx_op*v2v; dv2y=Dy_op*v2v; dw2x=Dx_op*w2v; dw2y=Dy_op*w2v;
            ibu2=1i*bb.*u2v; ibv2=1i*bb.*v2v; ibw2=1i*bb.*w2v;
            ibu3=1i*bc.*u3v; ibv3=1i*bc.*v3v; ibw3=1i*bc.*w3v;
            er3=EcG/Re_;
            cLT=er3*(lamT_b.*T1cv); cMT2=2*er3*(muT_b.*T1cv); cMT=er3*(muT_b.*T1cv);
            % a-slot d/dT_a: drop the T_a factor (forcing is coef*T_a*bracket)
            Gda = er3*lamT_b.*( du2x.*du3x + dv2y.*dv3y + ibw2.*ibw3 + 2*du2x.*dv3y + 2*du2x.*ibw3 + 2*dv2y.*ibw3 ) ...
                + 2*er3*muT_b.*( du2x.*du3x + dv2y.*dv3y + ibw2.*ibw3 ) ...
                + er3*muT_b.*( du2y.*du3y + 2*du2y.*dv3x + dv2x.*dv3x ) ...
                + er3*muT_b.*( dw2y.*dw3y + 2*dw2y.*ibv3 + ibv2.*ibv3 ) ...
                + er3*muT_b.*( dw2x.*dw3x + 2*dw2x.*ibu3 + ibu2.*ibu3 );
            Dbu = dg(cLT.*(du3x+2*dv3y+2*ibw3) + cMT2.*du3x)*Dx_op + dg(cMT.*(du3y+2*dv3x))*Dy_op + (1i*bb)*dg(cMT.*ibu3);
            Dbv = dg(cLT.*(dv3y+2*ibw3) + cMT2.*dv3y)*Dy_op + dg(cMT.*dv3x)*Dx_op + (1i*bb)*dg(cMT.*ibv3);
            Dbw = dg(cMT.*(dw3y+2*ibv3))*Dy_op + dg(cMT.*(dw3x+2*ibu3))*Dx_op + (1i*bb)*dg((cLT+cMT2).*ibw3);
            Dcu = dg((cLT+cMT2).*du2x)*Dx_op + dg(cMT.*du2y)*Dy_op + (1i*bc)*dg(cMT.*(2*dw2x+ibu2));
            Dcv = dg(cLT.*(dv2y+2*du2x) + cMT2.*dv2y)*Dy_op + dg(cMT.*(2*du2y+dv2x))*Dx_op + (1i*bc)*dg(cMT.*(2*dw2y+ibv2));
            Dcw = dg(cMT.*dw2y)*Dy_op + dg(cMT.*dw2x)*Dx_op + (1i*bc)*dg(cLT.*(ibw2+2*du2x+2*dv2y) + cMT2.*ibw2);
            if cba > 0
                Tcol_a=(cba-1)*bsize + T_col_local;
                [Ih,Jh_cols,Vh,Ia,Ja_cols,Va]=push_entries(Ih,Jh_cols,Vh,Ia,Ja_cols,Va, T_row, Tcol_a, Gda, conj_a);
            end
            if cbb > 0
                col_off_b=(cbb-1)*bsize;  bops={Dbu,0; Dbv,ny; Dbw,2*ny};
                for s=1:3
                    [oi,oj,ov]=find(bops{s,1}); if isempty(oi); continue; end
                    [Ih,Jh_cols,Vh,Ia,Ja_cols,Va]=push_entries(Ih,Jh_cols,Vh,Ia,Ja_cols,Va, ...
                        gmap(oi,gbase,row_off,4*ny), gmap(oj,gbase,col_off_b,bops{s,2}), ov, conj_b);
                end
            end
            if cbc > 0
                col_off_c=(cbc-1)*bsize;  cops={Dcu,0; Dcv,ny; Dcw,2*ny};
                for s=1:3
                    [oi,oj,ov]=find(cops{s,1}); if isempty(oi); continue; end
                    [Ih,Jh_cols,Vh,Ia,Ja_cols,Va]=push_entries(Ih,Jh_cols,Vh,Ia,Ja_cols,Va, ...
                        gmap(oi,gbase,row_off,4*ny), gmap(oj,gbase,col_off_c,cops{s,2}), ov, conj_c);
                end
            end

            % --- cubic pressure-work: Ec*ikM2 * u_b . grad(rho_a T_c) (E-row) ---
            ba=ctx.betavec(a); dr1x=Dx_op*r1; dr1y=Dy_op*r1; eik=EcG*ikM2g;
            if cba > 0
                col_off_a=(cba-1)*bsize;
                Pa = dg(eik*u2v.*T3v)*Dx_op + dg(eik*v2v.*T3v)*Dy_op ...
                   + dg(eik*(u2v.*dT3x + v2v.*dT3y + w2v.*(1i*(ba+bc)).*T3v));
                [oi,oj,ov]=find(Pa);
                if ~isempty(oi)
                    [Ih,Jh_cols,Vh,Ia,Ja_cols,Va]=push_entries(Ih,Jh_cols,Vh,Ia,Ja_cols,Va, ...
                        gmap(oi,gbase,row_off,4*ny), gmap(oj,gbase,col_off_a,3*ny), ov, conj_a);
                end
            end
            if cbb > 0
                col_off_b=(cbb-1)*bsize;
                [Ih,Jh_cols,Vh,Ia,Ja_cols,Va]=push_entries(Ih,Jh_cols,Vh,Ia,Ja_cols,Va, T_row, col_off_b+u_col_local, eik*(dr1x.*T3v + r1.*dT3x), conj_b);
                [Ih,Jh_cols,Vh,Ia,Ja_cols,Va]=push_entries(Ih,Jh_cols,Vh,Ia,Ja_cols,Va, T_row, col_off_b+v_col_local, eik*(dr1y.*T3v + r1.*dT3y), conj_b);
                [Ih,Jh_cols,Vh,Ia,Ja_cols,Va]=push_entries(Ih,Jh_cols,Vh,Ia,Ja_cols,Va, T_row, col_off_b+w_col_local, eik*(1i*ba.*r1.*T3v + r1.*(1i*bc.*T3v)), conj_b);
            end
            if cbc > 0
                col_off_c=(cbc-1)*bsize;
                Pc = dg(eik*u2v.*r1)*Dx_op + dg(eik*v2v.*r1)*Dy_op ...
                   + dg(eik*(u2v.*dr1x + v2v.*dr1y + w2v.*r1.*(1i*(ba+bc))));
                [oi,oj,ov]=find(Pc);
                if ~isempty(oi)
                    [Ih,Jh_cols,Vh,Ia,Ja_cols,Va]=push_entries(Ih,Jh_cols,Vh,Ia,Ja_cols,Va, ...
                        gmap(oi,gbase,row_off,4*ny), gmap(oj,gbase,col_off_c,4*ny), ov, conj_c);
                end
            end
          end
        end
    end
    % Fold this mode's triplets into the accumulators: bounds peak memory, and
    % resetting the cells per mode keeps each append O(1) rather than O(n^2).
    [Jh_acc,Ih,Jh_cols,Vh]=fold_tr(Jh_acc,Ih,Jh_cols,Vh,Ntot);
    [Ja_acc,Ia,Ja_cols,Va]=fold_tr(Ja_acc,Ia,Ja_cols,Va,Ntot);
end

%% Finalize: buffer-zone attenuation and zeroing the inflow rows
Jh = Jh_acc;
Ja = Ja_acc;

% Buffer-zone attenuation: scale each row by bufcf at its station.
ib = ctx.Bound.ibuff;
if ib <= nx
    % Build a diagonal scale matrix of size Ntot, with the per-station
    % bufcf value repeated across each station's rows.
    bufcf = ctx.Bound.bufcf(:);              % nx x 1
    per_station = repelem(bufcf, nv*ny);     % nv*ny per station
    scale_one_mode = per_station;            % bsize x 1
    scale_all = repmat(scale_one_mode, nphys, 1);
    D = spdiags(scale_all, 0, Ntot, Ntot);
    Jh = D * Jh;
    Ja = D * Ja;
end

% Zero the inflow station rows (i=1), as NLT_HNS leaves them: the gh-slot diagonal
% uses u/v/w_row, which includes station 1, so clear the whole first block.
inflow_rows_one = 1:nv*ny;   % all vars at station 1
for kk = 1:nphys
    r = (kk-1)*bsize + inflow_rows_one;
    Jh(r, :) = 0;
    Ja(r, :) = 0;
end

end


%% Local helpers

function [acc, I, J, V] = fold_tr(acc, I, J, V, Ntot)
%FOLD_TR  Concatenate the triplet cells (vertcat, O(n)) and fold them into the
%   sparse accumulator, which sums duplicate (row,col) entries, then clear.
if ~isempty(V)
    acc = acc + sparse(vertcat(I{:}), vertcat(J{:}), vertcat(V{:}), Ntot, Ntot);
    I = {}; J = {}; V = {};
end
end

function idx = gmap(gi, gbase, base_off, sub)
%GMAP  Map grid indices gi (1..ny*nx) to global matrix indices, using the
%   precomputed gbase lookup: idx = base_off + sub + gbase(gi). This is
%   equivalent to base_off + floor((gi-1)/ny)*nv*ny + sub + mod(gi-1,ny) + 1,
%   but O(1) per entry instead of recomputing floor/mod each time.
idx = base_off + sub + gbase(gi);
end

function [Ih, Jh, Vh, Ia, Ja, Va] = push_entries(Ih, Jh, Vh, Ia, Ja, Va, rows, cols, vals, to_antilin)
%PUSH_ENTRIES  Append (rows, cols, vals) as one cell to the holomorphic or
%   antilinear triplet-cell list. O(1) per call; the lists get concatenated
%   once per mode, in fold_tr.
rows = rows(:); cols = cols(:); vals = vals(:);
keep = vals ~= 0;
if ~any(keep); return; end
rows = rows(keep); cols = cols(keep); vals = vals(keep);
if to_antilin
    Ia{end+1} = rows;  Ja{end+1} = cols;  Va{end+1} = vals;
else
    Ih{end+1} = rows;  Jh{end+1} = cols;  Vh{end+1} = vals;
end
end

function [Dx_op, Dy_op] = build_diff_ops(ctx)
%BUILD_DIFF_OPS  Sparse single-mode physical derivative operators.
ny = ctx.ny; nx = ctx.nx;

% 1D xi derivative matrix (nx x nx): sparse Fornberg weights on the
% (possibly non-uniform) xi grid, standard sign convention. Reduces to a
% uniform 4th-order finite-difference (FD) stencil on a uniform grid.
Dxi_1D = ctx.D_xi;

% Kronecker product to lift this to the full ny*nx state (column-major,
% y fastest)
Dxi_op = kron(Dxi_1D, speye(ny));

% Wall-normal derivative (D1, size ny x ny) applied at each station:
% kron(I_nx, D1)
D1 = ctx.D1;
if ~issparse(D1); D1 = sparse(D1); end
Deta_op = kron(speye(nx), D1);

% Metric-weighted physical derivatives
xix_v  = ctx.xix(:);
etax_v = ctx.etax(:);
xiy_v  = ctx.xiy(:);
etay_v = ctx.etay(:);
Dx_op = spdiags(xix_v, 0, ny*nx, ny*nx) * Dxi_op + ...
        spdiags(etax_v,0, ny*nx, ny*nx) * Deta_op;
Dy_op = spdiags(xiy_v, 0, ny*nx, ny*nx) * Dxi_op + ...
        spdiags(etay_v,0, ny*nx, ny*nx) * Deta_op;
end

function [Dxx_op, Dyy_op, Dxy_op] = build_diff2_ops(ctx)
%BUILD_DIFF2_OPS  Sparse physical 2nd-derivative operators, matching phys_d2:
%   d_xx = xix^2 Xixi + 2 xix etax Xieta + etax^2 Etaeta + xixx Xi + etaxx Eta
%   d_yy and d_xy follow the same pattern, with their own metric coefficients.
ny = ctx.ny;  nx = ctx.nx;  n = ny*nx;
D_xi = ctx.D_xi;  D1 = ctx.D1;  if ~issparse(D1), D1 = sparse(D1); end
G = ctx.StabGrid;
D2 = G.D2;  if ~issparse(D2), D2 = sparse(D2); end
if isfield(G,'D_xi2') && ~isempty(G.D_xi2), D_xi2 = G.D_xi2; else, D_xi2 = D_xi*D_xi; end

Iny = speye(ny);  Inx = speye(nx);
Xixi   = kron(D_xi2, Iny);   % d^2/dxi^2
Etaeta = kron(Inx, D2);      % d^2/deta^2
Xieta  = kron(D_xi, D1);     % d/deta d/dxi
Xi     = kron(D_xi, Iny);    % d/dxi
Eta    = kron(Inx, D1);      % d/deta

dg = @(f) spdiags(f(:), 0, n, n);
Dxx_op = dg(G.xix.^2)*Xixi + dg(2*G.xix.*G.etax)*Xieta + dg(G.etax.^2)*Etaeta ...
       + dg(G.xixx)*Xi + dg(G.etaxx)*Eta;
Dyy_op = dg(G.xiy.^2)*Xixi + dg(2*G.xiy.*G.etay)*Xieta + dg(G.etay.^2)*Etaeta ...
       + dg(G.xiyy)*Xi + dg(G.etayy)*Eta;
Dxy_op = dg(G.xix.*G.xiy)*Xixi + dg(G.xix.*G.etay + G.xiy.*G.etax)*Xieta ...
       + dg(G.etax.*G.etay)*Etaeta + dg(G.xixy)*Xi + dg(G.etaxy)*Eta;
end
