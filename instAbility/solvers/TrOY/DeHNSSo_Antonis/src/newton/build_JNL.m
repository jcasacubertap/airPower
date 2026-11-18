function [Jh, Ja] = build_JNL(q_complex, ctx)
%BUILD_JNL  Analytical Jacobian of NLT_HNS w.r.t. q (2D case only).
%
%   [Jh, Ja] = build_JNL(q_complex, ctx)
%
%   Returns sparse complex matrices Jh, Ja of size (nphys*bsize) x (nphys*bsize)
%   representing the Wirtinger derivatives of the NLT forcing:
%
%       d f_NLT = Jh * dq_phys + Ja * conj(dq_phys)
%
%   Jh holds contributions where the differentiated slot is a physical mode.
%   Ja holds contributions where the differentiated slot is a conjugate of a
%   physical mode (hence the derivative acts on conj(q_phys)).
%
%   Only the momentum rows (u, v; w is zero in 2D) are filled. Inflow (i=1)
%   rows are zero. Buffer attenuation from ctx.Bound.bufcf is applied.
%
%   Assumptions: 2D TS case (beta = 0, w = 0 everywhere). The iu*beta*w_gh
%   terms in NLT_HNS vanish, so only the convective bilinear terms remain:
%
%       xm_j = - u_gh * Dx(u_jk) - v_gh * Dy(u_jk)
%       ym_j = - u_gh * Dx(v_jk) - v_gh * Dy(v_jk)
%
%   ctx fields used:
%     nphys, nf, ny, nx, bsize, phys_modes
%     HB, dxi, D1, xix, etax, xiy, etay, Bound.bufcf
%   Cached on first call (optional):
%     ctx.Dx_op, ctx.Dy_op  — one-mode physical derivative operators

nphys = ctx.nphys;
nf    = ctx.nf;
ny    = ctx.ny;
nx    = ctx.nx;
bsize = ctx.bsize;
phys  = ctx.phys_modes;
half  = round(nf/2);
Ntot  = nphys * bsize;

% One-mode physical derivative operators (ny*nx, column-major y-first)
if isfield(ctx, 'Dx_op') && ~isempty(ctx.Dx_op)
    Dx_op = ctx.Dx_op;
    Dy_op = ctx.Dy_op;
else
    [Dx_op, Dy_op] = build_diff_ops(ctx);
end

% Unpack all physical+conjugate modes into ny x nx arrays (u, v, w)
% w is included even in the 2D case because NLT_HNS still computes the
% w-momentum forcing zm_j and its derivatives are nonzero whenever w != 0.
u_mode = zeros(nf, ny, nx);
v_mode = zeros(nf, ny, nx);
w_mode = zeros(nf, ny, nx);
for k = 1:nphys
    j  = phys(k);
    jc = nf + 1 - j;
    off = (k-1)*bsize;
    qmat = reshape(q_complex(off+1:off+bsize), 4*ny, nx);
    u_mode(j,:,:) = qmat(1:ny,        :);
    v_mode(j,:,:) = qmat(ny+1:2*ny,   :);
    w_mode(j,:,:) = qmat(2*ny+1:3*ny, :);
    if jc ~= j
        u_mode(jc,:,:) = conj(qmat(1:ny,        :));
        v_mode(jc,:,:) = conj(qmat(ny+1:2*ny,   :));
        w_mode(jc,:,:) = conj(qmat(2*ny+1:3*ny, :));
    end
end

% Preallocate triplet accumulators
est_nnz_h = 0;
est_nnz_a = 0;
Ih = []; Jh_cols = []; Vh = [];
Ia = []; Ja_cols = []; Va = [];

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

% Mode-row insertion: adds entries at row block kk (physical idx in 1..nphys)
% to either Jh or Ja.

% Operator constants
Zx = sparse(ny*nx, ny*nx);

% Accumulate per physical mode
for kk = 1:nphys
    j = phys(kk);
    [gharray, jkarray] = find(ctx.HB(:,:,j));
    if isempty(gharray); continue; end

    % Row block offset for this mode's f_j
    row_off = (kk-1)*bsize;
    row_u   = row_off + (0:ny*nx-1);
    % Row index of u(y,i) within this mode's block = (i-1)*4*ny + y
    [Y, X] = ndgrid(1:ny, 1:nx);
    u_row_local = (X(:)-1)*4*ny +        Y(:);
    v_row_local = (X(:)-1)*4*ny +   ny + Y(:);
    w_row_local = (X(:)-1)*4*ny + 2*ny + Y(:);
    u_row = row_off + u_row_local;
    v_row = row_off + v_row_local;
    w_row = row_off + w_row_local;

    % Column-block row mapping: same layout
    u_col_local = u_row_local;
    v_col_local = v_row_local;
    w_col_local = w_row_local;

    for ip = 1:numel(gharray)
        gh = gharray(ip);
        jk = jkarray(ip);

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

        % ============================================================
        % Contributions from the gh slot (diagonal in y,x):
        %   ∂xm_j/∂u_gh = -diag(Dx_u_jk);   ∂xm_j/∂v_gh = -diag(Dy_u_jk)
        %   ∂ym_j/∂u_gh = -diag(Dx_v_jk);   ∂ym_j/∂v_gh = -diag(Dy_v_jk)
        %
        % These place diagonal entries at rows u_row/v_row and columns
        % u_col_gh / v_col_gh in column block "col_blk(gh)", target Jh or Ja.
        cb_gh  = col_blk(gh);
        conj_gh = is_conj(gh);
        if cb_gh > 0   % gh mode is active
            col_off_gh = (cb_gh-1)*bsize;
            ucol_gh = col_off_gh + u_col_local;
            vcol_gh = col_off_gh + v_col_local;
            wcol_gh = col_off_gh + w_col_local;

            [Ih, Jh_cols, Vh, Ia, Ja_cols, Va] = push_entries( ...
                Ih, Jh_cols, Vh, Ia, Ja_cols, Va, ...
                u_row, ucol_gh, -Dx_u_jk(:), conj_gh);
            [Ih, Jh_cols, Vh, Ia, Ja_cols, Va] = push_entries( ...
                Ih, Jh_cols, Vh, Ia, Ja_cols, Va, ...
                u_row, vcol_gh, -Dy_u_jk(:), conj_gh);
            [Ih, Jh_cols, Vh, Ia, Ja_cols, Va] = push_entries( ...
                Ih, Jh_cols, Vh, Ia, Ja_cols, Va, ...
                v_row, ucol_gh, -Dx_v_jk(:), conj_gh);
            [Ih, Jh_cols, Vh, Ia, Ja_cols, Va] = push_entries( ...
                Ih, Jh_cols, Vh, Ia, Ja_cols, Va, ...
                v_row, vcol_gh, -Dy_v_jk(:), conj_gh);
            [Ih, Jh_cols, Vh, Ia, Ja_cols, Va] = push_entries( ...
                Ih, Jh_cols, Vh, Ia, Ja_cols, Va, ...
                w_row, ucol_gh, -Dx_w_jk(:), conj_gh);
            [Ih, Jh_cols, Vh, Ia, Ja_cols, Va] = push_entries( ...
                Ih, Jh_cols, Vh, Ia, Ja_cols, Va, ...
                w_row, vcol_gh, -Dy_w_jk(:), conj_gh);
        end

        % ============================================================
        % Contributions from the jk slot (operator blocks):
        %   ∂xm_j/∂u_jk = -diag(u_gh)*Dx - diag(v_gh)*Dy
        %   ∂ym_j/∂v_jk = -diag(u_gh)*Dx - diag(v_gh)*Dy
        % xm has no v_jk; ym has no u_jk.
        cb_jk   = col_blk(jk);
        conj_jk = is_conj(jk);
        if cb_jk > 0   % jk mode is active
            col_off_jk = (cb_jk-1)*bsize;

            B_op = -spdiags(u_gh(:), 0, ny*nx, ny*nx) * Dx_op ...
                   -spdiags(v_gh(:), 0, ny*nx, ny*nx) * Dy_op;

            [bi, bj, bv] = find(B_op);
            if ~isempty(bi)
                rows_xm = row_off + (floor((bi-1)/ny))*4*ny + mod(bi-1,ny) + 1;
                cols_uj = col_off_jk + (floor((bj-1)/ny))*4*ny + mod(bj-1,ny) + 1;
                [Ih, Jh_cols, Vh, Ia, Ja_cols, Va] = push_entries( ...
                    Ih, Jh_cols, Vh, Ia, Ja_cols, Va, ...
                    rows_xm, cols_uj, bv, conj_jk);

                rows_ym = row_off + (floor((bi-1)/ny))*4*ny + ny + mod(bi-1,ny) + 1;
                cols_vj = col_off_jk + (floor((bj-1)/ny))*4*ny + ny + mod(bj-1,ny) + 1;
                [Ih, Jh_cols, Vh, Ia, Ja_cols, Va] = push_entries( ...
                    Ih, Jh_cols, Vh, Ia, Ja_cols, Va, ...
                    rows_ym, cols_vj, bv, conj_jk);

                rows_zm = row_off + (floor((bi-1)/ny))*4*ny + 2*ny + mod(bi-1,ny) + 1;
                cols_wj = col_off_jk + (floor((bj-1)/ny))*4*ny + 2*ny + mod(bj-1,ny) + 1;
                [Ih, Jh_cols, Vh, Ia, Ja_cols, Va] = push_entries( ...
                    Ih, Jh_cols, Vh, Ia, Ja_cols, Va, ...
                    rows_zm, cols_wj, bv, conj_jk);
            end
        end
    end
end

Jh = sparse(Ih, Jh_cols, Vh, Ntot, Ntot);
Ja = sparse(Ia, Ja_cols, Va, Ntot, Ntot);

% Buffer attenuation: scale each row by bufcf at its station.
ib = ctx.Bound.ibuff;
if ib <= nx
    % Build a diagonal scale of size Ntot with per-station bufcf.
    bufcf = ctx.Bound.bufcf(:);              % nx x 1
    per_station = repelem(bufcf, 4*ny);      % 4*ny per station
    scale_one_mode = per_station;            % bsize x 1
    scale_all = repmat(scale_one_mode, nphys, 1);
    D = spdiags(scale_all, 0, Ntot, Ntot);
    Jh = D * Jh;
    Ja = D * Ja;
end

% Inflow station rows (i=1) must be zero — mirrors NLT_HNS writing f(j,4*ny+1:end)
% and leaving f(j,1:4*ny) untouched (zero).
% Already handled implicitly: NLT never writes to station 1, and our row
% indices include station 1 contributions only via the operator rows
% (rows_xm / rows_ym) for x >= 2 cases. Actually, the "gh-slot diagonal" part
% uses u_row/v_row which INCLUDES station 1. Zero those rows.
inflow_rows_one = [ (1:ny), ny + (1:ny), 2*ny + (1:ny) ];   % u, v, w at station 1
for kk = 1:nphys
    r = (kk-1)*bsize + inflow_rows_one;
    Jh(r, :) = 0;
    Ja(r, :) = 0;
end

end


%% Local helpers ---------------------------------------------------------

function [Ih, Jh, Vh, Ia, Ja, Va] = push_entries(Ih, Jh, Vh, Ia, Ja, Va, rows, cols, vals, to_antilin)
%PUSH_ENTRIES  Append (rows, cols, vals) to holomorphic or antilinear list.
rows = rows(:); cols = cols(:); vals = vals(:);
keep = vals ~= 0;
if ~any(keep); return; end
rows = rows(keep); cols = cols(keep); vals = vals(keep);
if to_antilin
    Ia = [Ia; rows];  Ja = [Ja; cols];  Va = [Va; vals];
else
    Ih = [Ih; rows];  Jh = [Jh; cols];  Vh = [Vh; vals];
end
end

function [Dx_op, Dy_op] = build_diff_ops(ctx)
%BUILD_DIFF_OPS  Sparse single-mode physical derivative operators.
ny = ctx.ny; nx = ctx.nx;

% 1D xi derivative matrix (nx x nx) — sparse Fornberg weights on the
% (possibly non-uniform) ξ grid. Standard sign convention; reduces to
% uniform 4th-order FD on uniform grids.
Dxi_1D = ctx.D_xi;

% Kronecker to lift to ny*nx state (column-major y-first)
Dxi_op = kron(Dxi_1D, speye(ny));

% Wall-normal (D1 ny x ny) applied at each station: kron(I_nx, D1)
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

function D = fd1d4o_matrix(nx, d)
%FD1D4O_MATRIX  Sparse matrix form of FD1d4o (row-wise 4th-order FD).
% Recall FD1d4o flips the sign of the central coefficients (decreasing x).
D = sparse(nx, nx);
% Interior: central 4th order
for i = 3:nx-2
    D(i, i-2) = -1/12/d;
    D(i, i-1) =  2/3 /d;
    D(i, i+1) = -2/3 /d;
    D(i, i+2) =  1/12/d;
end
% Columns 1..2: backward (the "first two" of the flipped grid)
%   DD(1:2) = (-25/12*D(1:2) + 4*D(2:3) - 3*D(3:4) + 4/3*D(4:5) - 1/4*D(5:6)) / (-d)
for i = 1:2
    D(i, i  ) = (-25/12) / (-d);
    D(i, i+1) = ( 4    ) / (-d);
    D(i, i+2) = (-3    ) / (-d);
    D(i, i+3) = ( 4/3  ) / (-d);
    D(i, i+4) = (-1/4  ) / (-d);
end
% Columns end-1..end: forward
%   DD(end-1:end) = (-25/12*D(end-1:end) + 4*D(end-2:end-1) ... ) / d
for i = nx-1:nx
    D(i, i  ) = (-25/12) / d;
    D(i, i-1) = ( 4    ) / d;
    D(i, i-2) = (-3    ) / d;
    D(i, i-3) = ( 4/3  ) / d;
    D(i, i-4) = (-1/4  ) / d;
end
end
