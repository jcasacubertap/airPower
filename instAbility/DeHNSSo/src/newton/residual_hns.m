function F_real = residual_hns(q_real, ctx)
%RESIDUAL_HNS  Pure F(q) residual of the nonlinear HNS system (real arithmetic).
%
%   State vector is REAL of length 2*nphys*bsize:
%       q_real = [real(q_phys); imag(q_phys)]
%   where q_phys stacks the physical modes only. This keeps the JFNK
%   linearisation well-posed: conjugate modes are filled inside as
%   conj(q_phys), and NLT depends on both q_phys and conj(q_phys), so the
%   complex Jacobian is antilinear. Splitting into real/imag parts makes
%   the operator genuinely linear on q_real, which is what GMRES needs.

nphys = ctx.nphys;
nf    = ctx.nf;
ny    = ctx.ny;
nx    = ctx.nx;
bsize = ctx.bsize;
phys  = ctx.phys_modes;
N     = nphys * bsize;

% Real -> complex physical state
q_complex = q_real(1:N) + 1i * q_real(N+1:2*N);

% Unpack physical q into full (nf x ny x nx) u,v,w with conjugates
tmp.u = zeros(nf, ny, nx);
tmp.v = zeros(nf, ny, nx);
tmp.w = zeros(nf, ny, nx);
tmp.betavec = ctx.betavec;

for k = 1:nphys
    j   = phys(k);
    jc  = nf + 1 - j;
    off = (k-1)*bsize;
    qk_mat = reshape(q_complex(off+1:off+bsize), 4*ny, nx);
    uj = qk_mat(1:ny,        :);
    vj = qk_mat(ny+1:2*ny,   :);
    wj = qk_mat(2*ny+1:3*ny, :);
    tmp.u(j,:,:) = uj;  tmp.v(j,:,:) = vj;  tmp.w(j,:,:) = wj;
    if jc ~= j
        tmp.u(jc,:,:) = conj(uj);
        tmp.v(jc,:,:) = conj(vj);
        tmp.w(jc,:,:) = conj(wj);
    end
end

% NLT forcing (computed only for physical modes inside NLT_HNS)
f_all = NLT_HNS(ctx.StabGrid, 1:nf, tmp, ctx.D1, ctx.HB);

% Buffer attenuation on NLT forcing
ib = ctx.Bound.ibuff;
if ib <= nx
    bufcf_vec = repelem(ctx.Bound.bufcf(ib:nx), 4*ny);
    idx = (ib-1)*4*ny+1 : nx*4*ny;
    f_all(:, idx) = f_all(:, idx) .* bufcf_vec;
end

% Per-physical-mode residual (complex), with continuation factor AF on BC
AF = ctx.AF;
F_complex = zeros(N, 1);
for k = 1:nphys
    j   = phys(k);
    off = (k-1)*bsize;
    qk  = q_complex(off+1:off+bsize);
    fk  = f_all(j,:).';
    rBC = AF * ctx.rBC{k};
    fk(rBC ~= 0) = 0;
    F_complex(off+1:off+bsize) = ctx.ML{k}*qk - rBC - fk;
end

% Complex -> real
F_real = [real(F_complex); imag(F_complex)];
end
