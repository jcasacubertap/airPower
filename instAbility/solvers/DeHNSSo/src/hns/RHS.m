function [R] = RHS(Bound,Grid,Mode,Opt,StabRes,j)
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
% Construct the right-hand side vector for mode j.
%
%   R = RHS(Bound, Grid, Mode, Opt, StabRes, j)
%
%   Assembles the RHS from:
%     - Station 1: inflow shape function (phiIC), scaled by amplitude factor AF.
%     - Stations 2..nx: inhomogeneous wall/freestream boundary conditions (BCs)
%       from Bound.bcw/bct. Grid convention: index 1 = freestream (y=H), index
%       ny = wall (y=0). Freestream values injected only for IH_DR components
%       (Opt.bc_top).
%
%   Output: R [N x 1]: RHS vector for the global linear system ML * q = R.

%% Inputs
% Name    size         unit   explanation
% Bound   (1,1) struct [-]    boundary conditions and outflow-buffer coefficients
% Grid    (1,1) struct [-]    grid sizes, coordinates and inverse metric terms
% Mode    (1,1) struct [-]    mode bookkeeping: nf, active set RunJ and the interaction tables
% Opt     (1,1) struct [-]    solver options; see the option index in docs/DeHNSSo_guide.md
% StabRes (1,1) struct [-]    solution: shape functions, modal amplitudes, wavenumbers and iteration history
% j       (1,1)        [-]    mode counter of the mode the right-hand side is built for

%% Outputs
% Name size           unit   explanation
% R    (nvar*ny*nx,1) [-]    right-hand side: inflow shape at station 1, inhomogeneous wall and freestream conditions downstream

    ny = Grid.ny;
    R  = 0*StabRes.R;
    is_compr = isfield(Opt,'compressible') && strcmpi(Opt.compressible,'on');
    nvar     = 4 + double(is_compr);
    n_bc     = 3 + double(is_compr);

    % Station 1: inflow
    if ismember(j, Mode.RunJ)
        R(1:nvar*ny, 1) = Opt.AF * StabRes.phiIC(j,:);
    end

    % Stations 2..nx: inhomogeneous BCs (zero by default, i.e. homogeneous BCs)
    if j <= size(Bound.bcw,3) && (any(Bound.bcw(:,:,j),'all') || any(Bound.bct(:,:,j),'all'))
        i_vec = 2:Grid.nx;
        base  = nvar*ny*(i_vec - 1);
        % Bound.bcw/bct are the PHYSICAL peak amplitude 2|û|; seed the complex
        % amplitude û = value/2 so the reconstructed field (post-process 2|û|) equals
        % the value the user set (same convention as the inflow A0; DeHNSSo_guide §5).
        % Same row map the LHS uses, so each channel's value lands on the row
        % whose BC type was consulted. The freestream and wall maps differ in the
        % compressible case (see bc_rows).
        B = bc_rows(ny, is_compr);
        for k = 1:n_bc
            if strcmpi(Opt.bc_top{k}, 'IH_DR')
                R(base + B.top(k)) = 0.5 * Bound.bct(k, i_vec, j);   % freestream
            end
            R(base + B.wall(k)) = 0.5 * Bound.bcw(k, i_vec, j);      % wall
        end
    end
end

