function print_config(StabGrid, Stab, Opt)
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
% Print run configuration and memory estimates.
%
%   print_config(StabGrid, Stab, Opt)
%
%   Reports the grid, spectral content, inflow, compressibility, and solver
%   settings for the run, then estimates the memory required by the LHS
%   operator and its LU factorization (the LU decomposition used to solve the
%   linear system) and warns when that estimate is large.
%

%% Inputs
% Name     size unit   explanation
% StabGrid      [-]    stability grid (dimensions, coordinates, flow parameters)
% Stab          [-]    disturbance spectrum settings (M, N, IC, A0, ...)
% Opt           [-]    run options (buffer, lu_mode, parfor, ...)

%% Outputs
% None.

%% Configuration summary
is_compr = isfield(Opt,'compressible') && strcmpi(Opt.compressible,'on');
nvar     = 4 + double(is_compr);

fprintf('  --------------------------------------------\n');
fprintf('  Grid\n');
fprintf('    nx = %d,  ny = %d,  DOF = %d\n', StabGrid.nx, StabGrid.ny, nvar*StabGrid.ny*StabGrid.nx);
R_in  = sqrt(StabGrid.x(1)  * StabGrid.x(1));   % R = sqrt(x)*R0, R0 = sqrt(x_inflow)
R_out = sqrt(StabGrid.x(end) * StabGrid.x(1));
fprintf('    x = [%.2f, %.2f],  dxi = %.4f\n', StabGrid.x(1), StabGrid.x(end), StabGrid.dxi);
fprintf('    R = [%.2f, %.2f]   (sqrt(Re_x), local Blasius Reynolds)\n', R_in, R_out);
fprintf('  Spectral content\n');
n_phys = Stab.M * (2*Stab.N+1);
fprintf('    M = %d,  N = %d  (%d mode%s)\n', Stab.M, Stab.N, n_phys, plural(n_phys));
fprintf('    omega_0 = %.6f,  beta_0 = %.6f\n', Stab.omega_0, Stab.beta_0);
fprintf('  Inflow\n');
a0_tag = '';
if isfield(Opt,'linear') && strcmpi(Opt.linear,'on')
    a0_tag = ' [linear: shape scale-invariant; A0_fund sets the amplitude]';
end
a0_val = 0; if ~isempty(Stab.A0); a0_val = max(Stab.A0); end
fprintf('    IC: %s,  phase ref: %s,  A0 = %.4e%s\n', Stab.IC, Stab.phaseRef, a0_val, a0_tag);
if isfield(Stab,'bcw') && any(Stab.bcw(:))
    fprintf('    Wall actuation: ON\n');
end
if isfield(Stab,'bct') && any(Stab.bct(:))
    fprintf('    Freestream forcing: ON\n');
end
if isfield(Opt,'linear') && strcmpi(Opt.linear,'on')
    fprintf('    Mode: LINEAR (NLT inactive)\n');
else
    fprintf('    Mode: NONLINEAR\n');
end
if is_compr
    fprintf('  Compressible (u,v,w,rho,T)\n');
    fprintf('    Ma = %.3f,  Pr = %.3f,  gamma = %.2f,  S* = %.3f\n', ...
        StabGrid.Ma, StabGrid.Pr, StabGrid.gamma, StabGrid.Sstar);
    fprintf('    Wall thermal BC: %s\n', Opt.bc_wall_thermal);
end
fprintf('  Solver\n');
fprintf('    Buffer: %s from %d%%,  Conv: %.1e\n', Opt.buffer, Opt.xb, Opt.Conv);
% Channel names come from bc_rows so they cannot drift from the rows the operator
% actually constrains: (u,w,rho,T) compressible, (u,v,w) incompressible. There is
% no v channel at the compressible freestream.
B     = bc_rows(StabGrid.ny, is_compr);
pairs = cell(1, numel(B.top_names));
for k = 1:numel(B.top_names)
    pairs{k} = sprintf('%s=%s', B.top_names{k}, Opt.bc_top{k});
end
fprintf('    BC top: %s\n', strjoin(pairs, ', '));
fprintf('    LU mode: %s\n', Opt.lu_mode);
if isfield(Opt,'parfor') && strcmpi(Opt.parfor,'on')
    fprintf('    Parfor: on\n');
end

%% Memory estimates
nx    = StabGrid.nx;  ny = StabGrid.ny;
bsize = nvar * ny;
N_dof = bsize * nx;
D1_is_sparse = isfield(StabGrid,'D1') && issparse(StabGrid.D1) && (nnz(StabGrid.D1) < 0.2*ny^2);
if D1_is_sparse
    lhs_nnz_est = (20 + 2) * N_dof;
    eta_lbl = 'fd';
else
    lhs_nnz_est = (54 * (ny/40) + 2) * N_dof;
    eta_lbl = 'cheb';
end
lhs_mem_est = lhs_nnz_est * 20 / 1e9;

b_half      = 2 * bsize;
lu_mem_one  = 2 * N_dof * b_half * 20 * 0.75 / 1e9;
n_unique    = Stab.M + 1;
if Stab.N > 0
    n_unique = (Stab.M + 1) * (Stab.N + 1);
end
if ~isfield(Opt, 'lu_max_cache'); Opt.lu_max_cache = 3; end

if strcmpi(Opt.parfor,'on') && license('test','Distrib_Computing_Toolbox') && ~isempty(ver('parallel'))
    pool = gcp('nocreate');
    if isempty(pool); nw = 0; else; nw = min(pool.NumWorkers, n_unique); end
else
    nw = 0;
end
parfor_overhead = nw * n_unique * lu_mem_one;

fprintf('  Memory estimates (%s):\n', eta_lbl);
fprintf('    LHS:                   ~%.1f GB\n', lhs_mem_est);
fprintf('    LU per factorisation:  ~%.1f GB\n', lu_mem_one);
fprintf('    LU full   (cache %d):  ~%.1f GB\n', n_unique, n_unique * lu_mem_one);
fprintf('    LU auto   (cache %d):  ~%.1f GB\n', Opt.lu_max_cache, Opt.lu_max_cache * lu_mem_one);
fprintf('    LU single (cache 1):   ~%.1f GB\n', lu_mem_one);
fprintf('    Total peak (single):   ~%.1f GB\n', lu_mem_one + lhs_mem_est);
if nw > 0
    fprintf('    parfor workers (%d):   ~%.1f GB overhead\n', nw, parfor_overhead);
    fprintf('    Total peak (parfor):   ~%.1f GB\n', n_unique * lu_mem_one + lhs_mem_est + parfor_overhead);
end

lu_mem_total = lu_mem_one;
if strcmp(Opt.lu_mode, 'full'); lu_mem_total = n_unique * lu_mem_one; end
if strcmpi(Opt.parfor,'on'); lu_mem_total = lu_mem_total + parfor_overhead; end
if lu_mem_total > 20
    warning('Estimated LU memory: %.1f GB. Consider lu_mode = ''auto'' or ''none'' if RAM is limited.', lu_mem_total);
end
fprintf('  --------------------------------------------\n');

end

%% Local helpers

function s = plural(n)
%PLURAL  Return the plural suffix 's' unless n == 1.
if n == 1; s = ''; else; s = 's'; end
end
