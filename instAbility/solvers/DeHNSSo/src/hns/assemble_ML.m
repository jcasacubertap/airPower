function ML = assemble_ML(omega, beta, MLA1, MLA2, MLA3, MLA4, MLA5, MFD_diff)
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
% Per-mode harmonic operator for one (omega, beta) pair.
%
%   ML = ASSEMBLE_ML(omega, beta, MLA1, MLA2, MLA3, MLA4, MLA5, MFD_diff)
%
%   MLA1       omega- and beta-independent block
%   MLA2       linear in omega
%   MLA3       linear in beta
%   MLA4       quadratic in beta
%   MLA5       i*beta*d_x cross-viscous block (compressible only; zero sparse
%              for the incompressible operator)
%   MFD_diff   restores no-slip on the (0,0) mean-flow-distortion mode
%
%   Called by solve_modes to build the operator each mode is factorized from,
%   and by backward_error to rebuild the same operator for the residual check.

%% Inputs
% Name     size              unit   explanation
% omega    (1,1)             [-]    angular frequency of the mode
% beta     (1,1)             [-]    spanwise wavenumber of the mode
% MLA1     nvar*ny*nx sparse [-]    mode-invariant operator block, zeroth order in omega and beta
% MLA2     nvar*ny*nx sparse [-]    mode-invariant operator block, first order in omega
% MLA3     nvar*ny*nx sparse [-]    mode-invariant operator block, first order in beta
% MLA4     nvar*ny*nx sparse [-]    mode-invariant operator block, second order in beta
% MLA5     nvar*ny*nx sparse [-]    mode-invariant operator block, mixed omega-beta term
% MFD_diff nvar*ny*nx sparse [-]    difference operator for the mean-flow-distortion block

%% Outputs
% Name size              unit   explanation
% ML   nvar*ny*nx sparse [-]    harmonic operator for the given (omega, beta) pair

if omega == 0 && beta == 0
    ML = MLA1 + MFD_diff;
elseif omega == 0
    ML = MLA1 + MLA3*beta + MLA4*beta^2 + MLA5*beta;
elseif beta == 0
    ML = MLA1 + MLA2*omega;
else
    ML = MLA1 + MLA2*omega + MLA3*beta + MLA4*beta^2 + MLA5*beta;
end
end
