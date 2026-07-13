function ctx = factor_precond_full(ctx, q_real)
%FACTOR_PRECOND_FULL  Build and LU-factor the full mode-coupled preconditioner.
%
%   Assembles M = blkdiag(ML_j) - Jh_NLT(q), where Jh_NLT is the holomorphic
%   part of the analytical NLT Jacobian (build_JNL). The antiholomorphic part
%   Ja_NLT is dropped from the preconditioner — GMRES still sees the full
%   operator via JFNK, the preconditioner is just an approximation. This
%   keeps memory tractable (one N x N complex sparse matrix instead of the
%   2N x 2N real form).
%
%   Sets ctx.LU_full and ctx.use_full_precond = true. Subsequent calls to
%   apply_Minv use the global factorisation rather than the per-mode blocks.
%
%   Cost: build_JNL ~5-10 s + assembly ~1 s + sparse LU ~30-90 s.

N = ctx.nphys * ctx.bsize;
q_complex = q_real(1:N) + 1i * q_real(N+1:2*N);

% Free anything large we don't need before this point
if isfield(ctx, 'LU_full'); ctx.LU_full = []; end

t_jnl = tic;
[Jh_NLT, ~] = build_JNL(q_complex, ctx);
fprintf('   factor_full: build_JNL %.1fs (nnz=%.1fM)\n', toc(t_jnl), nnz(Jh_NLT)/1e6);

% Threshold tiny entries (< 1e-14 of max abs) — shaves ~20-30% of nnz
% without affecting the preconditioner's quality.
mx = max(abs(nonzeros(Jh_NLT)));
if mx > 0
    Jh_NLT = Jh_NLT .* (abs(Jh_NLT) > 1e-14*mx);
end

% Assemble M = blkdiag(ML_j) - Jh_NLT, freeing ML_block as soon as done
t_asm = tic;
ML_block = blkdiag(ctx.ML{:});
M_full = ML_block - Jh_NLT;
clear ML_block Jh_NLT
fprintf('   factor_full: assemble %.1fs (nnz(M)=%.1fM)\n', toc(t_asm), nnz(M_full)/1e6);

% Free the per-mode LU cache to recover several GB for the global factor.
ctx.LU = cell(ctx.nphys, 1);

% Incomplete LU as preconditioner — O(nnz(M)) memory, constructed in seconds.
% Try ILU(0) first (no fill); if it breaks down, fall back to a thresholded
% pivoted variant which is more robust at the cost of some extra fill.
t_lu = tic;
L = []; U = [];
try
    setup0 = struct('type','nofill');
    [L, U] = ilu(M_full, setup0);
    kind = 'ilu0';
catch ME0
    fprintf('   factor_full: ilu(0) failed (%s), trying ilutp ...\n', ME0.message);
end
if isempty(L)
    try
        setup1 = struct('type','ilutp','droptol',1e-4,'thresh',0.1);
        [L, U] = ilu(M_full, setup1);
        kind = 'ilutp';
    catch ME1
        error('factor_precond_full:iluFailed', ...
              'Both ILU(0) and ILUTP(1e-4) broke down: %s', ME1.message);
    end
end
fprintf('   factor_full: %s %.1fs (L+U: %.0f MB)\n', kind, toc(t_lu), (nnz(L)+nnz(U))*16/1e6);
clear M_full

% Store with identity row/col permutations so apply_Minv can use the same
% generic dispatch (it applies p on input and q on output).
N_total = ctx.nphys * ctx.bsize;
id = (1:N_total).';
ctx.LU_full = struct('L', L, 'U', U, 'p', id, 'q', id);
ctx.use_full_precond = true;
end
