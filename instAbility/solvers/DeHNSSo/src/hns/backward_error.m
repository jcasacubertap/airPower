function [eta_max, eta_mode, res_mode] = backward_error(Bound, StabGrid, Mode, Opt, Stab, StabRes, MLA1, MLA2, MLA3, MLA4, MLA5, MFD_diff)
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
% Componentwise backward error of the per-mode linear solves.
%
%   [omega_max, omega_mode, res_mode] = BACKWARD_ERROR(...)
%
%   Puts the solution each back-end returned back into the system it was handed
%   and reports how far off it is, mode by mode:
%
%       r   = ML*q - b
%       eta = max_i |r|_i / ( |ML|*|q| + |b| )_i     (Oettli-Prager, i not a BC row)
%       res = ||r||_2 / ||b||_2
%
%   eta is the componentwise backward error: the smallest relative perturbation
%   of each entry of ML and b for which q is the exact solution. The
%   componentwise form is the one that suits this operator, whose rows differ in
%   scale by many decades (second-derivative rows carry 1/dy^2), so a normwise
%   measure would be set by the largest row alone and would report a backward
%   error far below machine epsilon.
%
%   The maximum skips the Dirichlet rows. Those hold a single entry with b_i = 0
%   and are imposed rather than solved, so their ratio is identically 1 as soon
%   as the solution there is roundoff rather than an exact zero, which would pin
%   the maximum at 1 for every back-end.
%
%   eta answers what a factorization is responsible for: whether it solved the
%   system it was given. Unlike a difference against another back-end it needs no
%   reference run, so every back-end including the reference gets a number. A
%   backward-stable factorization returns eta at machine epsilon whatever the
%   conditioning of ML; a large eta means the factorization itself is at fault.
%
%   Call immediately after solve_modes, before AVAL, Anderson or update_modes
%   touch StabRes.q or Stab.f: at that point q is the raw solve output and the
%   right-hand side can be rebuilt exactly as the solve saw it. RHS reads only
%   StabRes.phiIC and Opt.AF, neither of which has changed yet.
%
%   eta_mode and res_mode are indexed by the physical modes, in the order
%   Mode.RunJ gives them; eta_max is the largest eta over those modes.

%% Inputs
% Name     size              unit   explanation
% Bound    (1,1) struct      [-]    boundary conditions and outflow-buffer coefficients
% StabGrid (1,1) struct      [-]    stability grid, coordinate metrics and base-flow fields on it
% Mode     (1,1) struct      [-]    mode bookkeeping: nf, active set RunJ and the interaction tables
% Opt      (1,1) struct      [-]    solver options; see the option index in docs/DeHNSSo_guide.md
% Stab     (1,1) struct      [-]    disturbance specification: modal truncation, frequencies, inflow condition and amplitudes
% StabRes  (1,1) struct      [-]    solution: shape functions, modal amplitudes, wavenumbers and iteration history
% MLA1     nvar*ny*nx sparse [-]    mode-invariant operator block, zeroth order in omega and beta
% MLA2     nvar*ny*nx sparse [-]    mode-invariant operator block, first order in omega
% MLA3     nvar*ny*nx sparse [-]    mode-invariant operator block, first order in beta
% MLA4     nvar*ny*nx sparse [-]    mode-invariant operator block, second order in beta
% MLA5     nvar*ny*nx sparse [-]    mode-invariant operator block, mixed omega-beta term
% MFD_diff nvar*ny*nx sparse [-]    difference operator for the mean-flow-distortion block

%% Outputs
% Name     size      unit   explanation
% eta_max  (1,1)     [-]    largest componentwise backward error over the physical modes
% eta_mode (1,nphys) [-]    componentwise backward error of each physical mode
% res_mode (1,nphys) [-]    relative residual norm ||r||/||b|| of each physical mode

phys_modes = Mode.RunJ(Mode.RunJ >= round(Mode.nf/2));
nphys      = numel(phys_modes);
eta_mode   = zeros(1, nphys);
res_mode   = zeros(1, nphys);

% Same Dirichlet masking of the nonlinear forcing that solve_modes applies, so
% the right-hand side matches the one that was actually solved.
Stab_f     = Stab.f;
bc_dir     = full(sum(MLA1 ~= 0, 2)) == 1;
bc_dir_mfd = full(sum((MLA1 + MFD_diff) ~= 0, 2)) == 1;
for jj = 1:size(Stab_f, 1)
    if Stab.Mvec(jj) == 0 && Stab.Nvec(jj) == 0
        Stab_f(jj, bc_dir_mfd) = 0;
    else
        Stab_f(jj, bc_dir) = 0;
    end
end

for k = 1:nphys
    j  = phys_modes(k);
    ML = assemble_ML(StabRes.omegavec(j), StabRes.betavec(j), ...
                     MLA1, MLA2, MLA3, MLA4, MLA5, MFD_diff);

    R_j = RHS(Bound, StabGrid, Mode, Opt, StabRes, j);
    f_j = Stab_f(j,:);
    f_j(R_j ~= 0) = 0;
    b = R_j + f_j.';

    x = StabRes.q(j,:).';
    r = ML*x - b;

    if Stab.Mvec(j) == 0 && Stab.Nvec(j) == 0
        keep = ~bc_dir_mfd;
    else
        keep = ~bc_dir;
    end
    den  = abs(ML)*abs(x) + abs(b);
    keep = keep & den > 0;
    eta_mode(k) = max(abs(r(keep)) ./ den(keep));
    res_mode(k) = norm(r) / max(norm(b), realmin);
end
eta_max = max(eta_mode);
end
