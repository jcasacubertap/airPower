function [StabRes, info] = solve_newton(Stab, StabGrid, Opt)
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
% Pseudo-arclength continuation in the inflow amplitude A0_fund.
%
%   [StabRes, info] = solve_newton(Stab, StabGrid, Opt)
%
%   Continues the nonlinear HNS solution along A0_fund = lambda from
%   Opt.newton_start to Opt.newton_target, for cases that will not converge from
%   a cold start at the target. Keller continuation: the predictor is a tangent
%   plus a curvature built from the two preceding tangents; the corrector is a
%   modified Newton on the coupled real-form Jacobian, factorised with MUMPS and
%   reused across steps until it ages out or stalls, a factorisation costing far
%   more than a back-solve. lambda enters F only through the inflow BC, hence
%   affinely,
%       F(q,lambda) = newton_residual(q,ctx0) + (lambda-lam0)*F_lambda,
%   with F_lambda finite-differenced once.
%
%   Opt (defaults): newton_target (required), newton_start (Stab.A0_fund),
%   newton_seed2 (1.04*start), newton_reuse (8), newton_it_thresh (8),
%   newton_maxit (25), newton_dsmax (3), newton_maxstep (400), newton_tol (1e-6),
%   newton_seedfile, newton_ckfile, newton_resume ('off').
%
%   See also NEWTON_CONTEXT, NLT_JACOBIAN, NEWTON_RESIDUAL.

%% Inputs
% Name     size         unit   explanation
% Stab     (1,1) struct [-]    disturbance specification: modal truncation, frequencies, inflow condition and amplitudes
% StabGrid (1,1) struct [-]    stability grid, coordinate metrics and base-flow fields on it
% Opt      (1,1) struct [-]    solver options; see the option index in docs/DeHNSSo_guide.md

%% Outputs
% Name    size         unit   explanation
% StabRes (1,1) struct [-]    solution: shape functions, modal amplitudes, wavenumbers and iteration history
% info    (1,1) struct [-]    continuation report: status, number of factorizations, step count and the branch trajectory

%% Options and continuation parameters
t_all = tic; cum = 0;

if ~isfield(Opt,'newton_target') || isempty(Opt.newton_target)
    error('solve_newton:noTarget','Opt.newton_target (target A0_fund) is required.');
end
A_target = Opt.newton_target;
lam0 = Stab.A0_fund; if isfield(Opt,'newton_start') && ~isempty(Opt.newton_start); lam0 = Opt.newton_start; end
lam1 = 1.04*lam0;    if isfield(Opt,'newton_seed2') && ~isempty(Opt.newton_seed2); lam1 = Opt.newton_seed2; end
% Defaults are the modified continuation of the paper; corrector='full' with
% predictor='tangent' recovers Dijkstra et al. (2014), Sec. 2.2.
corrector = lower(get_o(Opt,'newton_corrector','modified'));
predictor = lower(get_o(Opt,'newton_predictor','curvature'));
if ~ismember(corrector, {'modified','full'})
    error('solve_newton:corrector', ...
          'Opt.newton_corrector must be ''modified'' or ''full'' (got ''%s'').', corrector);
end
if ~ismember(predictor, {'curvature','tangent'})
    error('solve_newton:predictor', ...
          'Opt.newton_predictor must be ''curvature'' or ''tangent'' (got ''%s'').', predictor);
end
full_newton = strcmp(corrector,'full');
use_curv    = strcmp(predictor,'curvature');
if full_newton && (isfield(Opt,'newton_reuse') || isfield(Opt,'newton_it_thresh'))
    warning('solve_newton:reuseIgnored', ...
            ['newton_reuse / newton_it_thresh are ignored when ' ...
             'newton_corrector = ''full'' (the factor is never reused).']);
end

Kmax        = get_o(Opt,'newton_reuse', 8);
it_thresh   = get_o(Opt,'newton_it_thresh', 8);
reuse_maxit = get_o(Opt,'newton_maxit', 25);
dsmax       = get_o(Opt,'newton_dsmax', 3);
maxstep     = get_o(Opt,'newton_maxstep', 400);
arctol      = get_o(Opt,'newton_tol', 1e-6);
% The seeds only finite-difference an initial tangent, which the corrector then
% repairs, so they need no arclength tolerance. 1e-4 saves a factorisation per
% seed; affects only the first run on a grid, since newton_seedfile covers reruns.
seedtol     = get_o(Opt,'newton_seed_tol', 1e-6);
% Arclength metric. 'absolute' (default) is ds^2 = |dq|^2 + dlam^2, ill-conditioned
% when the state norm dwarfs the parameter, so ds tracks the state and lambda
% crawls. 'relative' is ds^2 = |dq|^2/|qc|^2 + (dlam/lc)^2.
arcmode = lower(get_o(Opt,'newton_arclength','absolute'));
is_rel  = strcmp(arcmode,'relative');
seedfile    = get_o(Opt,'newton_seedfile', []);
ckfile      = get_o(Opt,'newton_ckfile', []);
do_resume   = isfield(Opt,'newton_resume') && strcmpi(Opt.newton_resume,'on');

fprintf('\n--- Pseudo-arclength continuation: A0 %.4e -> %.4e (u_rms %.4f) ---\n', ...
        lam0, A_target, A_target/sqrt(2));

% Base flow fields for newton_context (as assembled by the callers/drivers).
BF.Re  = StabGrid.Re;
BF.U   = StabGrid.U;   BF.V   = StabGrid.V;   BF.W   = StabGrid.W;
BF.dxU = StabGrid.dxU; BF.dxV = StabGrid.dxV; BF.dxW = StabGrid.dxW;
BF.dyU = StabGrid.dyU; BF.dyV = StabGrid.dyV; BF.dyW = StabGrid.dyW;

% LHS_compressible reads BF.Ma, rho, T and their derivatives, as in main.m.
if isfield(Opt,'compressible') && strcmpi(Opt.compressible,'on')
    req = {'rho','T','dxrho','dyrho','dxT','dyT','dxxU','dxxV','dxxW','dxxT', ...
           'ddyU','ddyV','ddyW','ddyT','dxyU','dxyV', ...
           'Ma','Pr','gamma','Sstar','mu','dmuT','d2muT2'};
    for kf = 1:numel(req); BF.(req{kf}) = StabGrid.(req{kf}); end
end

%% Context at a given lambda
% Nested function: any name shared with the parent IS the parent's variable, so
% these are kept distinct. Reusing qc/q0/Mode would overwrite the continuation state.
    function [ctx_o, qv_o, SG_o, Bnd_o, Mode_o, Opt_o, Stab_o, SR_o] = mkctx(lam)
        S_l = Stab; S_l.A0_fund = lam;
        [SG_o, Init_l, Mode_o, Opt_o, Stab_o, SR_o] = init(StabGrid, Opt, S_l);
        [Mode_o, SR_o] = inflow(Mode_o, Stab_o, SG_o, SR_o, Opt_o);
        [Bnd_o]        = setup_boundary(SG_o, Opt_o, Stab_o);
        Mode_o.RunJ    = 1:Mode_o.nf; Mode_o.Aactive = ones(1, Mode_o.nf);
        ctx_o = newton_context(BF, Bnd_o, SG_o, Init_l, Mode_o, Opt_o, Stab_o, SR_o);
        ctx_o.AF = 1;
        Nl_l = ctx_o.nphys*ctx_o.bsize; qc_l = zeros(Nl_l,1);
        for k_l = 1:ctx_o.nphys
            j_l = ctx_o.phys_modes(k_l); off_l = (k_l-1)*ctx_o.bsize;
            qc_l(off_l+1:off_l+ctx_o.bsize) = SR_o.q(j_l,:).';
        end
        qv_o = [real(qc_l); imag(qc_l)];
    end

fprintf(' building ctx0 @ %.4e ...\n', lam0);
[ctx0, q_init, SG0, Bnd0, Mode0, Opt0, Stab0, SR0] = mkctx(lam0);
N = ctx0.nphys*ctx0.bsize; twoN = 2*N;
ML = blkdiag(ctx0.ML{:});
fprintf(' N(complex)=%d  2N=%d  nphys=%d bsize=%d nf=%d\n', ...
        N, twoN, ctx0.nphys, ctx0.bsize, ctx0.nf);

% The MUMPS ordering must match the wall-normal discretisation or factorisation
% fails with INFOG(1)=-13 on the other type: AMD for Chebyshev (dense D1), METIS
% for FD (sparse D1). newton_mumps_ordering overrides.
if ~isfield(Opt,'newton_mumps_ordering') || isempty(Opt.newton_mumps_ordering)
    is_fd = isfield(ctx0,'D1') && issparse(ctx0.D1);
    Opt.newton_mumps_ordering = tern(is_fd, 7, 0);   % 7=auto(->METIS on FD), 0=AMD
    fprintf(' MUMPS ordering: %s (ICNTL(7)=%d), auto-selected from %s wall-normal operator\n', ...
            tern(is_fd,'auto/METIS','AMD'), Opt.newton_mumps_ordering, tern(is_fd,'sparse FD','dense Cheb'));
end

%% dF/dlambda (affine in lambda: inflow BC only)
h = 1e-6; [ctxh, ~] = mkctx(lam0 + h);
Fa = newton_residual(q_init, ctx0); Fb = newton_residual(q_init, ctxh);
% phi = -g, with g = dF/dlambda of the paper. Carrying the sign here is what
% turns the paper's dx = z1 - dA0*z2 into the addition dq = a + dl*b below.
Flam = (Fb - Fa)/h; phi = -Flam; phi = phi(:); Flam = Flam(:);
[ctx2, ~] = mkctx(lam0 + 50*h); Fchk = newton_residual(q_init, ctx2);
fprintf(' ||F_lambda||=%.3e  linearity rel.err=%.2e\n', ...
        norm(Flam), norm(Fchk - (Fa + 50*h*Flam))/max(norm(Fchk), eps));
clear ctxh ctx2

evalF = @(q,lam) newton_residual(q, ctx0) + (lam - lam0)*Flam;

    function Jr = Jreal(q)
        qcx = q(1:N) + 1i*q(N+1:twoN);
        [Jh, Ja] = nlt_jacobian(qcx, ctx0); A = ML - Jh; B = -Ja;
        Jr = [real(A)+real(B), -imag(A)+imag(B); imag(A)+imag(B), real(A)-real(B)];
    end
% Natural-parameter Newton for the seeds. Locals are suffixed _n: nested, so a
% bare `id` would be the parent's MUMPS handle, destroyed here every iteration.
    function [qs_, ls_] = arcscale(qc_, lc_, rel_)
    % Arclength scales. 'absolute' returns 1,1; 'relative' returns the current
    % state and parameter magnitudes, so ds measures relative change in each.
        if rel_
            qs_ = norm(qc_); if ~(qs_ > 0), qs_ = 1; end
            ls_ = abs(lc_);  if ~(ls_ > 0), ls_ = 1; end
        else
            qs_ = 1; ls_ = 1;
        end
    end

    function q_out = natnewton(q_in, lam, tol, mx)
        q_out = q_in;
        for it_n = 1:mx
            F_n = evalF(q_out, lam); nF_n = norm(F_n);
            fprintf('   nat it%d lam=%.5e ||F||=%.3e\n', it_n, lam, nF_n);
            if nF_n < tol, break; end
            Jr_n = Jreal(q_out); id_n = mumps_solve('factor', Jr_n, Opt, 0, lam);
            [dq_n, id_n] = mumps_solve('bsolve', id_n, Jr_n, -F_n, Opt);
            id_n = mumps_solve('destroy', id_n);
            q_out = q_out + dq_n;
        end
    end

%% Two seed points
% Skipped when resuming: the resume branch rebuilds the tangent from the
% checkpoint, and computing seeds anyway means two cold-start solves at the
% resume amplitude, the regime this method exists to avoid.
have_ck = do_resume && ~isempty(ckfile) && exist(ckfile,'file') == 2;
if have_ck
    fprintf('\n -- seed points skipped (resuming from %s) --\n', ckfile);
    q0 = q_init; q1 = q_init; tq = zeros(size(q_init)); tl = 1; ds = 0;
    qc = q_init(:); lc = lam0;
else
    fprintf('\n -- seed points --\n'); t_seed = tic;
    if ~isempty(seedfile) && exist(seedfile,'file') == 2
        L = load(seedfile); q0 = L.q0(:); q1 = L.q1(:);
        fprintf(' loaded seed from %s\n', seedfile);
    else
        q0 = natnewton(q_init, lam0, seedtol, 30);
        q1 = natnewton(q0,     lam1, seedtol, 30);
        if ~isempty(seedfile)
            save(seedfile, 'q0', 'q1', 'lam0', 'lam1', '-v7.3');
            fprintf(' saved seed to %s\n', seedfile);
        end
    end
    fprintf(' [timing] seed: %.1f min\n', toc(t_seed)/60);

    dq = q1 - q0; dl = lam1 - lam0;
    qc = q1(:); lc = lam1;
    [qs, ls] = arcscale(qc, lc, is_rel);
    nrm = sqrt((dq'*dq)/qs^2 + (dl/ls)^2);
    tq = (dq/qs)/nrm; tl = (dl/ls)/nrm; ds = max(nrm, 3e-5); tq = tq(:);
    fprintf(' initial tangent: dA0/ds=%.3e  ds=%.3e\n', tl, ds);
end

%% Cross-step reuse state
% age = 1e9 forces a factorisation on the first step.
id = []; Jr = []; age = 1e9; it_prev = 0; nfac_total = 0;
traj = []; tq_old = []; tl_old = []; ds_prev = 0.35;

if do_resume && ~isempty(ckfile) && exist(ckfile,'file') == 2
    L = load(ckfile); qc = L.q(:); lc = L.l; traj = L.traj;
    if isfield(L,'nfac_total'); nfac_total = L.nfac_total; end
    fprintf('   [factorising: resume] ... '); tfac = tic;
    Jr = Jreal(qc); id = mumps_solve('factor', Jr, Opt, size(traj,1), lc); nfac_total = nfac_total + 1; age = 0;
    fprintf('done in %.1f min\n', toc(tfac)/60);
    [w, id] = mumps_solve('bsolve', id, Jr, phi, Opt);
    [qs, ls] = arcscale(qc, lc, is_rel);
    w = w*ls/qs;                       % dq/dlam -> d(q/qs)/d(lam/ls)
    % Next tangent, the paper's [-z2', 1]' normalised, signed so lambda increases.
    tl = 1/sqrt(w'*w + 1); tq = tl*w; tq = tq(:);
    if tl < 0, tl = -tl; tq = -tq; end
    ds = get_o(Opt,'newton_ds', 2.0);
    fprintf(' RESUMED from %s: A0=%.6e u_rms=%.5e, %d steps done, dA0/ds=%.3e, ds=%.3e\n', ...
            ckfile, lc, lc/sqrt(2), size(traj,1), tl, ds);
end

info.status = 'maxstep'; info.nfac_total = 0; info.steps = 0;

%% Main continuation loop
for step = (size(traj,1)+1):maxstep
    t_step = tic;
    [qs, ls] = arcscale(qc, lc, is_rel);
    if tl > 0 && lc + ds*ls*tl > A_target, ds = (A_target - lc)/(ls*tl); end
    if use_curv && ~isempty(tq_old)
        cq = (tq - tq_old)/ds_prev; cl = (tl - tl_old)/ds_prev;
        qp = qc + qs*(ds*tq + 0.5*ds^2*cq); lp = lc + ls*(ds*tl + 0.5*ds^2*cl);
    else
        qp = qc + qs*ds*tq; lp = lc + ls*ds*tl;
    end
    q = qp(:); l = lp; conv = false; nfac_step = 0;

    % A non-finite predictor must reject the step, not reach mumps_factor and throw
    % Not a step-length problem: the usual cause is 0.5*ds^2*cq with a rotating tangent.
    if ~all(isfinite(q)) || ~isfinite(l)
        if use_curv && ~isempty(tq_old)
            fprintf(' STEP %d: non-finite predictor; retrying first-order (ds=%.3e)\n', step, ds);
            q = qc + ds*tq; q = q(:); l = lc + ds*tl;
        end
        if ~all(isfinite(q)) || ~isfinite(l)
            fprintf(' STEP %d: predictor still non-finite; ds %.3e -> %.3e\n', step, ds, ds/2);
            if ~isempty(id), id = mumps_solve('destroy', id); id = []; end
            ds = ds/2;
            if ds < 1e-6
                fprintf(' ds too small; stop\n'); info.status = 'ds_underflow'; break
            end
            tq_old = []; tl_old = [];     % drop the curvature history that caused it
            continue
        end
    end

    if full_newton
        % Dijkstra Eq. (2.12) evaluates the Jacobian at the current iterate, so it is
        % refactorised every corrector iteration and z2 re-solved. b stays empty if
        % the predictor already converged; the tangent then carries over.
        b = [];
        for it = 1:reuse_maxit
            q = q(:); F = evalF(q, l); F = F(:);
            Nc  = tq(:).'*((q - qc)/qs) + tl*((l - lc)/ls) - ds;
            res = sqrt(F'*F + Nc^2);
            if res < arctol, conv = true; break; end
            fprintf('   [factorising: full Newton, iteration %d] ... ', it); tfac = tic;
            if ~isempty(id), id = mumps_solve('destroy', id); id = []; end
            Jr = Jreal(q); id = mumps_solve('factor', Jr, Opt, step, l);
            nfac_step = nfac_step + 1; nfac_total = nfac_total + 1; age = 0;
            fprintf('done in %.1f min\n', toc(tfac)/60);
            [a, id] = mumps_solve('bsolve', id, Jr, -F, Opt);  a = a(:);
            [b, id] = mumps_solve('bsolve', id, Jr, phi, Opt); b = b(:);
            dl = (-Nc - tq(:).'*a/qs)/(tq(:).'*b/qs + tl/ls); dl = dl(1);
            dq = a + dl*b; q = q(:) + dq(:); l = l + dl;
        end
    else
    % reuse the standing factor, or refactor at the predictor
    need_fac = isempty(id) || age >= Kmax || it_prev > it_thresh;
    if need_fac
        if isempty(id)
            why = 'no factor';
        elseif age >= Kmax
            why = sprintf('age %d >= reuse %d', age, Kmax);
        else
            why = sprintf('last corrector took %d its > %d', it_prev, it_thresh);
        end
        fprintf('   [factorising: %s] ... ', why); tfac = tic;
        if ~isempty(id), id = mumps_solve('destroy', id); id = []; end
        Jr = Jreal(q); id = mumps_solve('factor', Jr, Opt, step, l);
        nfac_step = 1; nfac_total = nfac_total + 1; age = 0;
        fprintf('done in %.1f min\n', toc(tfac)/60);
    end
    % b is the paper's z2 (up to the sign in phi), the solve against the parameter
    % derivative. It is fixed within the corrector, so it is taken once per step.
    [b, id] = mumps_solve('bsolve', id, Jr, phi, Opt); b = b(:);

    % Corrector, Eq. (newton_two_solves) of the paper: Nc is the constraint S, a is
    % z1 and b is z2, and each iteration combines the two so the step satisfies S=0.
    for it = 1:reuse_maxit
        q = q(:); F = evalF(q, l); F = F(:);
        Nc  = tq(:).'*((q - qc)/qs) + tl*((l - lc)/ls) - ds;
        res = sqrt(F'*F + Nc^2);
        if res < arctol, conv = true; break; end
        [a, id] = mumps_solve('bsolve', id, Jr, -F, Opt); a = a(:);
        dl = (-Nc - tq(:).'*a/qs)/(tq(:).'*b/qs + tl/ls); dl = dl(1);
        dq = a + dl*b; q = q(:) + dq(:); l = l + dl;
    end

    % stale-factor fallback: refactor at current q, finish with a fresh factor
    if ~conv
        fprintf('   [factorising: corrector stalled after %d its] ... ', reuse_maxit); tfac = tic;
        if ~isempty(id), id = mumps_solve('destroy', id); id = []; end
        Jr = Jreal(q); id = mumps_solve('factor', Jr, Opt, step, l);
        nfac_step = nfac_step + 1; nfac_total = nfac_total + 1; age = 0;
        fprintf('done in %.1f min\n', toc(tfac)/60);
        [b, id] = mumps_solve('bsolve', id, Jr, phi, Opt); b = b(:);
        for it = 1:40
            q = q(:); F = evalF(q, l); F = F(:);
            Nc  = tq(:).'*((q - qc)/qs) + tl*((l - lc)/ls) - ds;
            res = sqrt(F'*F + Nc^2);
            if res < arctol, conv = true; break; end
            [a, id] = mumps_solve('bsolve', id, Jr, -F, Opt); a = a(:);
            dl = (-Nc - tq(:).'*a/qs)/(tq(:).'*b/qs + tl/ls); dl = dl(1);
            dq = a + dl*b; q = q(:) + dq(:); l = l + dl;
        end
    end
    end

    dt = toc(t_step)/60; cum = cum + dt;
    urms = l/sqrt(2); amp = norm(q(1:N) + 1i*q(N+1:twoN));
    fprintf(' STEP %d: A0=%.6e u_rms=%.5e |q|=%.4e it=%d fac=%d age=%d %s (ds=%.3e) [%.1f min, cum %.1f, nfac_tot=%d]\n', ...
            step, l, urms, amp, it, nfac_step, age, tern(conv,'ok','FAIL'), ds, dt, cum, nfac_total);

    if ~conv
        if ~isempty(id), id = mumps_solve('destroy', id); id = []; end
        ds = ds/2;
        if ds < 1e-6, fprintf(' ds too small; stop\n'); info.status = 'ds_underflow'; break; end
        continue
    end
    qc = q; lc = l; it_prev = it; age = age + 1;
    traj(end+1,:) = [l, urms, amp, it, nfac_step];
    if ~isempty(ckfile)
        save(ckfile, 'q', 'l', 'N', 'twoN', 'step', 'traj', 'nfac_total', '-v7.3');
    end

    if l >= A_target - 1e-9
        % LAND ON THE TARGET. The predictor is clipped to A_target, but pseudo-arclength
        % lets the corrector move lambda too, so the accepted state normally sits past
        % it, and the overshoot goes straight into the answer: the fundamental responds
        % ~linearly to the inflow amplitude and the harmonics ~quadratically. The
        % overshot state is a good guess, so this takes 1-2 iterations, but natnewton
        % refactorises each one.
        if l > A_target*(1 + 1e-6)
            fprintf('\n landing on target: A0 %.6e -> %.6e (overshoot %.2f%%)\n', ...
                    l, A_target, 100*(l/A_target - 1));
            % Try a frozen-Jacobian landing on the factor already held: fixed lambda needs
            % no bordering, so each iteration is one back-solve. natnewton stays the
            % fallback, since the held factor may be Kmax steps stale and the landing moves
            % lambda by the whole overshoot.
            land_ok = false;
            if ~full_newton && ~isempty(id)
                q_try = q(:);
                for it_l = 1:reuse_maxit
                    F_l = evalF(q_try, A_target); F_l = F_l(:);
                    if norm(F_l) < arctol, land_ok = true; break; end
                    [dq_l, id] = mumps_solve('bsolve', id, Jr, -F_l, Opt);
                    q_try = q_try + dq_l(:);
                end
                if land_ok
                    q = q_try;
                    fprintf('   landed on held factor in %d back-solves (no refactorisation)\n', it_l);
                end
            end
            if ~land_ok
                fprintf('   held factor too stale; falling back to natnewton\n');
                q = natnewton(q, A_target, arctol, reuse_maxit);
            end
            l = A_target;
            % StabRes is unpacked from qc (the last ACCEPTED step) and
            % info.lambda from lc, not from q/l, so the landed state must be
            % promoted here or it is computed and then silently discarded.
            qc = q; lc = l;
            amp = norm(q(1:N) + 1i*q(N+1:twoN));
            traj(end+1,:) = [l, l/sqrt(2), amp, NaN, 0];
            if ~isempty(ckfile)
                save(ckfile, 'q', 'l', 'N', 'twoN', 'step', 'traj', 'nfac_total', '-v7.3');
            end
        end
        if ~isempty(id), id = mumps_solve('destroy', id); id = []; end
        fprintf('\n >>> reached target A0=%.6e u_rms=%.5f |q|=%.4e\n', l, l/sqrt(2), amp);
        fprintf(' >>> TOTAL wall time: %.1f min (%.2f h) over %d steps | nfac_total=%d\n', ...
                toc(t_all)/60, toc(t_all)/3600, step, nfac_total);
        info.status = 'converged';
        break
    end

    if ~isempty(b)
        tq_old = tq; tl_old = tl; ds_prev = ds;
        w = b*ls/qs; tln = 1/sqrt(w'*w + 1); tqn = tln*w;
        if (tqn'*tq + tln*tl) < 0, tqn = -tqn; tln = -tln; end
        tq = tqn; tl = tln;
    end
    % Step-size control. The original rule grew ds at it<=6 and cut it only at
    % it>=reuse_maxit-2 (=23), leaving a dead band over it=7..22 where a strained
    % corrector produced no response, letting the next predictor go non-finite.
    % iterations, ds was left alone, and the next predictor went non-finite. React
    % inside the band instead, growing only on easy steps and cutting early on strain.
    it_easy = get_o(Opt,'newton_it_easy', 4);
    it_hard = get_o(Opt,'newton_it_hard', 10);
    if     it <= it_easy,          ds = min(ds*1.2, dsmax);
    elseif it >= reuse_maxit-2,    ds = ds*0.6;
    elseif it >= it_hard,          ds = ds*0.8;     % strain: back off before failing
    end
end
if ~isempty(id), id = mumps_solve('destroy', id); end

info.nfac_total = nfac_total;
info.steps      = size(traj,1);
info.traj       = traj;
info.wall_min   = toc(t_all)/60;
info.lambda     = lc;

%% Unpack final state into StabRes
StabRes = SR0; Mode = Mode0;
q_complex = qc(1:N) + 1i*qc(N+1:twoN);
for k = 1:ctx0.nphys
    j   = ctx0.phys_modes(k);
    off = (k-1)*ctx0.bsize;
    qj  = q_complex(off+1:off+ctx0.bsize).';
    StabRes.q(j, :) = qj;
    jc = ctx0.nf + 1 - j;
    if jc ~= j; StabRes.q(jc, :) = conj(qj); end
end
[StabRes, ~, ~] = update_modes(SG0, Bnd0, Mode, Opt0, Stab0, StabRes);

% Apply the SAME postprocessing as the Picard path (main.m) so both solvers
% return StabRes in one convention: shape functions normalised to unit max|u'|,
% amplitudes A = 2|u'|_max (MFD = |u'|_max), conjugate modes dropped. Uses the
% internal, fully-active Mode0; main.m skips its own postprocess on this path.
Init_pp.iu = 1i;
StabRes = postprocess(SG0, Init_pp, Mode, Opt0, StabRes);

fprintf('\n=== solve_newton DONE. %.1f min (%.2f h). nfac_total=%d, steps=%d, status=%s ===\n', ...
        info.wall_min, info.wall_min/60, nfac_total, info.steps, info.status);
end

%% Local helpers
function v = get_o(Opt, name, dflt)
if isfield(Opt, name) && ~isempty(Opt.(name)); v = Opt.(name); else; v = dflt; end
end
function s = tern(c, a, b)
if c, s = a; else, s = b; end
end
