function ctx = newton_context(BF, Bound, StabGrid, Init, Mode, Opt, Stab, StabRes)
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
% Assemble the Newton solver context.
%
%   ctx = newton_context(BF, Bound, StabGrid, Init, Mode, Opt, Stab, StabRes)
%
%   Builds the per-mode LHS operators, their block-preconditioner LU factors and
%   the inhomogeneous boundary-condition vectors, packed into one context struct.
%   Used by solve_newton (pseudo-arclength continuation in A0_fund) and by the
%   nlt_jacobian instrumentation in paper/checks/.

%% Inputs
% Name     size         unit   explanation
% BF       (1,1) struct [-]    base flow on its own grid
% Bound    (1,1) struct [-]    boundary conditions and outflow-buffer coefficients
% StabGrid (1,1) struct [-]    stability grid, coordinate metrics and base-flow fields on it
% Init     (1,1) struct [-]    inflow initialization data
% Mode     (1,1) struct [-]    mode bookkeeping: nf, active set RunJ and the interaction tables
% Opt      (1,1) struct [-]    solver options; see the option index in docs/DeHNSSo_guide.md
% Stab     (1,1) struct [-]    disturbance specification: modal truncation, frequencies, inflow condition and amplitudes
% StabRes  (1,1) struct [-]    solution: shape functions, modal amplitudes, wavenumbers and iteration history

%% Outputs
% Name size         unit   explanation
% ctx  (1,1) struct [-]    Newton context: operator blocks, harmonic-balancing registry, mode bookkeeping and sizes

%% Assemble LHS operators
is_compr = isfield(Opt,'compressible') && strcmpi(Opt.compressible,'on');
nvar     = 4 + is_compr;   % 4 incompressible | 5 compressible
fprintf(' building LHS (%s) ...', ternary_c(is_compr,'compressible','incompressible'));
t = tic;
MLA5 = [];
if is_compr
    [MLAMFD, MLA1, MLA2, MLA3, MLA4, MLA5] = LHS_compressible(BF, Bound, StabGrid, Init, Stab, StabRes, Opt);
else
    [MLAMFD, MLA1, MLA2, MLA3, MLA4]       = LHS(BF, Bound, StabGrid, Init, Stab, StabRes, Opt);
end
% Helmholtz-resonator reactive wall (Opt.hr), written into MLA1 only, so the (0,0)
% MFD mode keeps no-slip through MFD_diff.
if isfield(Bound,'hr') && ~isempty(Bound.hr)
    MLA1 = apply_resonator_wall(MLA1, Bound.hr, BF, StabGrid, is_compr);
end
% Keep only the difference: the MFD (0,0) operator differs from MLA1 in its
% boundary rows alone, so it is rebuilt below as MLA1 + MFD_diff.
MFD_diff = MLAMFD - MLA1;
clear MLAMFD

fprintf(' %.1fs\n', toc(t));

%% Active modes and context fields
% Conjugate modes are stored in the lower half of the spectrum, so the physical
% ones are those at or above the midpoint.
half  = round(Mode.nf/2);
phys  = Mode.RunJ(Mode.RunJ >= half);   % active physical modes only
nphys = numel(phys);

% Sizes. bsize is the length of one mode's state, so the full Newton state is
% nphys*bsize complex, twice that in the real form the solver works in.
ctx.phys_modes = phys;
ctx.nphys      = nphys;
ctx.nf         = Mode.nf;
ctx.ny         = StabGrid.ny;
ctx.nx         = StabGrid.nx;
ctx.nvar       = nvar;
ctx.bsize      = nvar*ctx.ny*ctx.nx;

% Grid and wall-normal operator, read by the NLT and by nlt_jacobian.
ctx.StabGrid   = StabGrid;
ctx.D1         = StabGrid.D1;

% HB is the mode-interaction table: HB(:,:,j) lists the (gh,jk) pairs whose
% product forces mode j, which is what nlt_jacobian walks.
ctx.HB         = Stab.HB;
ctx.omegavec   = StabRes.omegavec;
ctx.betavec    = StabRes.betavec;

% Inputs carried through untouched, so the residual can be re-evaluated at any
% lambda without the caller passing them again.
ctx.Bound      = Bound;
ctx.Mode       = Mode;
ctx.BF_orig    = BF;
ctx.Init       = Init;
ctx.Stab       = Stab;
ctx.StabRes    = StabRes;
ctx.Opt        = Opt;

% Spectral indices, and the curvilinear metrics used to form physical-space
% derivatives from the (xi,eta) operators.
ctx.Mvec       = Stab.Mvec;
ctx.Nvec       = Stab.Nvec;
ctx.dxi        = StabGrid.x(2) - StabGrid.x(1);   % legacy scalar, for display only
ctx.D_xi       = StabGrid.D_xi;                    % sparse non-uniform FD matrix in xi
ctx.xix        = StabGrid.xix;
ctx.etax       = StabGrid.etax;
ctx.xiy        = StabGrid.xiy;
ctx.etay       = StabGrid.etay;

% Per-mode operator and boundary-condition right-hand side, filled below.
ctx.ML  = cell(nphys, 1);
ctx.rBC = cell(nphys, 1);

Opt_rhs = Opt;
Opt_rhs.AF = 1;   % store UNIT amplitude BC; newton_residual scales by ctx.AF
key_cache = struct();

%% Per-mode LHS blocks
% The LHS is polynomial in (omega, beta), so each mode's operator is combined from
% the precomputed blocks rather than reassembled. Modes sharing an (omega, beta)
% pair share one operator through key_cache.
for k = 1:nphys
    j     = phys(k);
    omega = StabRes.omegavec(j);
    beta  = StabRes.betavec(j);
    key   = ml_key(omega, beta);

    if ~isfield(key_cache, key)
        % MLA5 is the i*beta*d_x cross-viscous block; compressible only, empty otherwise.
        b5 = 0; if is_compr && ~isempty(MLA5); b5 = beta; end
        % The branches are the one polynomial with its zero terms dropped, so a 2D
        % or steady mode does not pay for sparse additions that vanish anyway.
        if omega == 0 && beta == 0
            ML = MLA1 + MFD_diff;              % (0,0) MFD: its own boundary rows
        elseif omega == 0
            ML = MLA1 + MLA3*beta + MLA4*beta^2;
            if b5, ML = ML + MLA5*b5; end
        elseif beta == 0
            ML = MLA1 + MLA2*omega;
        else
            ML = MLA1 + MLA2*omega + MLA3*beta + MLA4*beta^2;
            if b5, ML = ML + MLA5*b5; end
        end
        key_cache.(key) = struct('ML', ML);
    end
    ctx.ML{k} = key_cache.(key).ML;
    % Inflow BC at unit amplitude; newton_residual scales it by the continuation AF.
    ctx.rBC{k} = RHS(Bound, StabGrid, Mode, Opt_rhs, StabRes, j);
end

end


function key = ml_key(omega, beta)
% Cache key for one (omega, beta) pair; the (0,0) mean-flow distortion is named
% separately since it carries its own boundary rows.
if omega == 0 && beta == 0
    key = 'MFD';
else
    key = matlab.lang.makeValidName(sprintf('o%g_b%g', omega, beta));
end
end

function s = ternary_c(c,a,b); if c, s=a; else, s=b; end; end
