function [Opt,StabRes] = AVAL(Grid,Opt,Mode,StabRes)
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
% Check the peak disturbance amplitude and damp the inflow if it is too large.
%
%   [Opt, StabRes] = AVAL(Grid, Opt, Mode, StabRes)
%
%   Finds the peak velocity amplitude of the fundamental mode and stores it in
%   Opt.AVAL. If that peak exceeds the ceiling Opt.AMAX, the amplitude factor
%   Opt.AF is rescaled so the peak comes back down to Opt.AMAX, and the whole
%   solution vector is scaled by the new Opt.AF.

%% Inputs
% Name    size         unit   explanation
% Grid    (1,1) struct [-]    grid sizes, coordinates and inverse metric terms
% Opt     (1,1) struct [-]    solver options; see the option index in docs/DeHNSSo_guide.md
% Mode    (1,1) struct [-]    mode bookkeeping: nf, active set RunJ and the interaction tables
% StabRes (1,1) struct [-]    solution: shape functions, modal amplitudes, wavenumbers and iteration history

%% Outputs
% Name    size         unit   explanation
% Opt     (1,1) struct [-]    solver options; see the option index in docs/DeHNSSo_guide.md
% StabRes (1,1) struct [-]    solution: shape functions, modal amplitudes, wavenumbers and iteration history

%% Reshape the solution vector into the field array phi
    ny_   = Grid.ny;
    nx_   = Grid.nx;
    nvar  = size(StabRes.phi, 2) / ny_;
    n_vel = 3;                              % u, v, w slots
    for j = Mode.RunJ
        StabRes.phi(j,:,:) = reshape(StabRes.q(j,:), nvar*ny_, nx_);
    end

%% Peak velocity amplitude of the fundamental mode
    [Opt.AVAL,~] = max(max(abs(squeeze(StabRes.phi(round(Mode.nf/2)+1,1:n_vel*ny_,:)))));

%% Damp the inflow amplitude if the peak exceeds the ceiling AMAX
    if Opt.AVAL >= Opt.AMAX
        Opt.AF = Opt.AMAX/Opt.AVAL;
    end

    StabRes.q = Opt.AF*StabRes.q;
end
