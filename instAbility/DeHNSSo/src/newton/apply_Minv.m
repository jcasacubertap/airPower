function dq_real = apply_Minv(r_real, ctx)
%APPLY_MINV  Preconditioner apply in real arithmetic.
%
%   Interprets r_real as [real(r_phys); imag(r_phys)], applies a complex
%   preconditioner, and returns [real(dq_phys); imag(dq_phys)].
%
%   Two preconditioner modes are supported:
%     ctx.use_full_precond = true  -> ctx.LU_full (single global N x N solve,
%       captures cross-mode coupling via Jh_NLT, see factor_precond_full.m)
%     otherwise -> per-mode block-diagonal solve via ctx.LU{k} (faster but
%       drops all NLT coupling, fine when nonlinearity is mild)

nphys = ctx.nphys;
bsize = ctx.bsize;
N     = nphys * bsize;

r_complex  = r_real(1:N) + 1i * r_real(N+1:2*N);
dq_complex = zeros(N, 1);

if isfield(ctx, 'use_full_precond') && ctx.use_full_precond
    s = ctx.LU_full;
    dq_complex(s.q) = s.U \ (s.L \ r_complex(s.p));
else
    for k = 1:nphys
        off  = (k-1)*bsize;
        r_k  = r_complex(off+1:off+bsize);
        lu_s = ctx.LU{k};
        x_k  = zeros(size(r_k));
        x_k(lu_s.q) = lu_s.U \ (lu_s.L \ r_k(lu_s.p));
        dq_complex(off+1:off+bsize) = x_k;
    end
end

dq_real = [real(dq_complex); imag(dq_complex)];
end
