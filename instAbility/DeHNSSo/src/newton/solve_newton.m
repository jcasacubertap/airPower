function [StabRes, Mode, info] = solve_newton(BF, Bound, StabGrid, Init, Mode, Opt, Stab, StabRes)
%SOLVE_NEWTON  Newton solver for the nonlinear HNS system.
%
%   JFNK Newton: matrix-free Krylov solve with per-mode block-LU
%   preconditioner (or full ILU via Opt.newton_precond='full').
%
%   When called as fallback from the convergence loop in main.m, uses the current
%   state (StabRes.q) as warm start and respects Mode.RunJ (active modes).
%
%   Opt fields (defaults):
%     newton_tol     1e-10    absolute ||F|| tolerance
%     newton_maxit   20       max Newton iterations
%     newton_precond 'block'  'block' = per-mode sparse LU preconditioner
%                             'full'  = single ILU of blkdiag(ML) - Jh_NLT
%                                       (lower memory, captures cross-mode NLT)

if ~isfield(Opt, 'newton_tol');     Opt.newton_tol     = 1e-10;   end
if ~isfield(Opt, 'newton_maxit');   Opt.newton_maxit   = 20;      end
if ~isfield(Opt, 'newton_precond'); Opt.newton_precond = 'block'; end
% Opt.newton_precond: 'block' — per-mode block LU (~nphys x per-LU, fastest GMRES)
%                     'full'  — ILU(0) of blkdiag(ML) - Jh_NLT (~2-10x less memory,
%                                captures cross-mode coupling; GMRES ~2x slower)

% When called directly (not from Picard fallback), activate all modes.
% When called from the convergence loop, Mode.RunJ is already set there.
half = round(Mode.nf/2);
n_phys_active = numel(Mode.RunJ(Mode.RunJ >= half));
if n_phys_active <= 1
    % Only fundamental active — we're called directly, not from Picard
    Mode.RunJ    = 1:Mode.nf;
    Mode.Aactive = ones(1, Mode.nf);
end

% Build LHS, per-mode LU, BC vectors for active modes only
ctx = build_ctx(BF, Bound, StabGrid, Init, Mode, Opt, Stab, StabRes);
N = ctx.nphys * ctx.bsize;

% Warm start from current StabRes (not linear solution — we're correcting Picard)
ctx.AF = Opt.AF;
q_complex = zeros(N, 1);
for k = 1:ctx.nphys
    j   = ctx.phys_modes(k);
    off = (k-1)*ctx.bsize;
    q_complex(off+1:off+ctx.bsize) = StabRes.q(j, :).';
end
q = [real(q_complex); imag(q_complex)];

% Newton iteration
F   = residual_hns(q, ctx);
nrm = norm(F);

fprintf(' init: ||F|| = %.3e\n', nrm);
info.history = nrm;
info.status  = 'maxit';
info.iters   = 0;

[q, nrm, info] = run_jfnk(q, F, nrm, ctx, Opt, N, info);

% Unpack
q_complex = q(1:N) + 1i * q(N+1:2*N);
for k = 1:ctx.nphys
    j   = ctx.phys_modes(k);
    off = (k-1)*ctx.bsize;
    qj  = q_complex(off+1:off+ctx.bsize).';
    StabRes.q(j, :) = qj;
    jc = ctx.nf + 1 - j;
    if jc ~= j; StabRes.q(jc, :) = conj(qj); end
end
[StabRes, Mode, ~] = update_modes(StabGrid, Bound, Mode, Opt, Stab, StabRes);
end


%% JFNK — matrix-free Newton with per-mode block-LU preconditioner
function [q, nrm, info] = run_jfnk(q, F, nrm, ctx, Opt, N, info)

fprintf(' %-4s %-11s %-11s %-5s %-8s %-6s\n', 'it','||F||','alpha','gm','t_J','t(s)');
ew_eta = 0.3;
krylov_maxit = 200;

% Build the full ILU preconditioner if requested — replaces the per-mode
% block LUs from build_ctx with a single ILU of M = blkdiag(ML) - Jh_NLT.
% Much lower memory (no fill), captures cross-mode coupling via Jh_NLT.
if isfield(Opt,'newton_precond') && strcmpi(Opt.newton_precond,'full')
    fprintf(' Building full ILU preconditioner ...\n');
    ctx = factor_precond_full(ctx, q);
end

% Per-mode ||F_k|| history for convergence plot
bsize = ctx.bsize;
nphys = ctx.nphys;
info.mode_F = zeros(nphys, 1);  % current per-mode ||F||
info.mode_F_history = [];       % nphys x n_iters
for kk = 1:nphys
    off = (kk-1)*bsize;
    info.mode_F(kk) = norm(F(off+1:off+bsize)) + norm(F(N+off+1:N+off+bsize));
end
info.mode_F_history = info.mode_F;

for k = 1:Opt.newton_maxit
    if nrm < Opt.newton_tol
        info.status = 'converged'; info.iters = k-1;
        fprintf(' converged at iter %d (||F|| = %.2e)\n', k-1, nrm);
        return
    end
    nrm_prev = nrm;
    t_it = tic;

    % Build analytical Jacobian (Jh, Ja) for fast matvec
    q_complex = q(1:N) + 1i * q(N+1:2*N);
    t_j = tic;
    [Jh, Ja] = build_JNL(q_complex, ctx);
    A_jac = blkdiag(ctx.ML{:}) - Jh;
    B_jac = -Ja;
    t_build = toc(t_j);

    Jop = @(v_real) jac_matvec(v_real, A_jac, B_jac, N);
    Mop = @(r) apply_Minv(r, ctx);

    % Use GMRES without restart (full Krylov subspace)
    [dq, ~, ~, gi] = gmres(Jop, -F, krylov_maxit, ew_eta, 1, Mop);
    total_gi = gi(2);

    % Armijo
    alpha = 1; accepted = false;
    for ls = 1:10
        q_trial = q + alpha * dq;
        F_trial = residual_hns(q_trial, ctx);
        nrm_trial = norm(F_trial);
        if isfinite(nrm_trial) && nrm_trial < (1 - 1e-4*alpha) * nrm_prev
            accepted = true; break
        end
        alpha = alpha * 0.5;
    end

    if accepted
        q = q_trial; F = F_trial; nrm = nrm_trial;
        ew_new = 0.9 * (nrm / nrm_prev)^2;
        sg = 0.9 * ew_eta^2;
        if sg > 0.1; ew_new = max(ew_new, sg); end
        ew_eta = min(0.5, max(1e-3, ew_new));
    end

    % Per-mode ||F_k||
    for kk = 1:nphys
        off = (kk-1)*bsize;
        info.mode_F(kk) = norm(F(off+1:off+bsize)) + norm(F(N+off+1:N+off+bsize));
    end
    info.mode_F_history(:, end+1) = info.mode_F;

    status = '.';
    if ~accepted; status = 'X'; end
    fprintf(' %-4d %-11.3e %-11.2g %-5d %-8.1f %-6.1f %s\n', ...
            k, nrm, alpha, total_gi, t_build, toc(t_it), status);
    info.history(end+1) = nrm;

    if ~isfinite(nrm); info.status = 'diverged'; info.iters = k; return; end

    % Early-bail if ‖F‖ is not materially decreasing — catches:
    % (a) rejected-step cycles (J near-singular → α → 0)
    % (b) plateau at discretization limit (‖F‖ stuck at some non-zero floor)
    % Stall = relative reduction < 1% over last 5 iters.
    n_window = 5;
    if k >= n_window
        recent = info.history(end-n_window+1:end);
        rel_drop = (max(recent) - min(recent)) / max(recent);
        if rel_drop < 1e-2
            fprintf([' Newton stalled: ||F|| = %.3e, ' ...
                     '%.2f%% drop over %d iters — aborting.\n'], ...
                     nrm, 100*rel_drop, n_window);
            info.status = 'stalled'; info.iters = k; return;
        end
    end
end
info.iters = Opt.newton_maxit;
end


%% jac_matvec — analytical J·v in real arithmetic
function w = jac_matvec(v_real, A, B, N)
%   A = ML - Jh (holomorphic), B = -Ja (antiholomorphic)
%   J_real * [vr;vi] = [Re(A*vc + B*conj(vc)); Im(A*vc + B*conj(vc))]
%   where vc = vr + i*vi
vc = v_real(1:N) + 1i * v_real(N+1:2*N);
wc = A * vc + B * conj(vc);
w  = [real(wc); imag(wc)];
end


%% build_ctx — LHS, per-mode LU, BC vectors (active modes only)
function ctx = build_ctx(BF, Bound, StabGrid, Init, Mode, Opt, Stab, StabRes)

fprintf(' building LHS ...');
t = tic;
[MLAMFD, MLA1, MLA2, MLA3, MLA4] = LHS(BF, Bound, StabGrid, Init, Stab, StabRes, Opt);
MFD_diff = MLAMFD - MLA1;
clear MLAMFD
fprintf(' %.1fs\n', toc(t));

half  = round(Mode.nf/2);
phys  = Mode.RunJ(Mode.RunJ >= half);   % active physical modes only
nphys = numel(phys);

ctx.phys_modes = phys;
ctx.nphys      = nphys;
ctx.nf         = Mode.nf;
ctx.ny         = StabGrid.ny;
ctx.nx         = StabGrid.nx;
ctx.bsize      = 4*ctx.ny*ctx.nx;
ctx.StabGrid   = StabGrid;
ctx.D1         = StabGrid.D1;
ctx.HB         = Stab.HB;
ctx.omegavec   = StabRes.omegavec;
ctx.betavec    = StabRes.betavec;
ctx.Bound      = Bound;
ctx.Mode       = Mode;
ctx.BF_orig    = BF;
ctx.Init       = Init;
ctx.Stab       = Stab;
ctx.StabRes    = StabRes;
ctx.Opt        = Opt;
ctx.Mvec       = Stab.Mvec;
ctx.Nvec       = Stab.Nvec;
ctx.dxi        = StabGrid.x(2) - StabGrid.x(1);   % legacy scalar (for display only)
ctx.D_xi       = StabGrid.D_xi;                    % sparse non-uniform FD matrix
ctx.xix        = StabGrid.xix;
ctx.etax       = StabGrid.etax;
ctx.xiy        = StabGrid.xiy;
ctx.etay       = StabGrid.etay;

ctx.ML  = cell(nphys, 1);
ctx.LU  = cell(nphys, 1);
ctx.rBC = cell(nphys, 1);

Opt_rhs = Opt;
Opt_rhs.AF = 1;   % store UNIT amplitude BC; residual_hns scales by ctx.AF
key_cache = struct();

% Skip per-mode LU factorisation when the full ILU preconditioner will be
% built later (factor_precond_full); it frees the per-mode LUs anyway.
% Skip per-mode block LU when newton_precond='full' (ILU preconditioner
% is built later in run_jfnk).
skip_block_lu = isfield(Opt,'newton_precond') && strcmpi(Opt.newton_precond,'full');

for k = 1:nphys
    j     = phys(k);
    omega = StabRes.omegavec(j);
    beta  = StabRes.betavec(j);
    key   = ml_key(omega, beta);

    if ~isfield(key_cache, key)
        if omega == 0 && beta == 0
            ML = MLA1 + MFD_diff;
        elseif omega == 0
            ML = MLA1 + MLA3*beta + MLA4*beta^2;
        elseif beta == 0
            ML = MLA1 + MLA2*omega;
        else
            ML = MLA1 + MLA2*omega + MLA3*beta + MLA4*beta^2;
        end
        if skip_block_lu
            fprintf('   ML mode (%d,%d) built (no block-LU, full ILU will be used)\n', ...
                    Stab.Mvec(j), Stab.Nvec(j));
            key_cache.(key) = struct('ML', ML, 'L', [], 'U', [], 'p', [], 'q', []);
        else
            fprintf('   LU mode (%d,%d) ...', Stab.Mvec(j), Stab.Nvec(j));
            t_lu = tic;
            [L, U, p, q_p] = lu(ML, 'vector');
            fprintf(' %.1fs (%.0f MB)\n', toc(t_lu), (nnz(L)+nnz(U))*16/1e6);
            key_cache.(key) = struct('ML', ML, 'L', L, 'U', U, 'p', p, 'q', q_p);
        end
    end
    entry = key_cache.(key);
    ctx.ML{k} = entry.ML;
    if ~skip_block_lu
        ctx.LU{k} = struct('L', entry.L, 'U', entry.U, 'p', entry.p, 'q', entry.q);
    end
    ctx.rBC{k} = RHS(Bound, StabGrid, Mode, Opt_rhs, StabRes, j);
end
end


function key = ml_key(omega, beta)
if omega == 0 && beta == 0
    key = 'MFD';
else
    key = matlab.lang.makeValidName(sprintf('o%g_b%g', omega, beta));
end
end
