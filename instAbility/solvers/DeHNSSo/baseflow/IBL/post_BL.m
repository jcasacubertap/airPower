function [BL] = post_BL(BL)
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
% Boundary-layer integral quantities.
%
%   BL = post_BL(BL)
%
%   Adds displacement/momentum thickness, shape factor H12, the 99% thickness,
%   and (compressible) the Stanton number to the BL struct.

%% Inputs
% Name size         unit   explanation
% BL   (1,1) struct [-]    boundary-layer profiles and derivatives

%% Outputs
% Name size         unit   explanation
% BL   (1,1) struct [-]    boundary-layer profiles and derivatives

fprintf('\nPost-processing BL solution\n')

[ny,nx]  = size(BL.u); % Matrix dimensions
[~,imax] = max(BL.y);  % Index of freestream

% Dependent properties
Uinf        = BL.Ue(1);
Winf        = BL.We;
BL.q        = sqrt(Uinf^2 + Winf^2);
BL.Cp       = 1-(BL.Ue/Uinf).^2;
BL.l        = sqrt(BL.S*BL.nu/Uinf);
BL.Re_l     = BL.l*Uinf/BL.nu;
BL.chord    = max(BL.x);

% Displacement & Momentum thickness
if BL.model
    rho  = BL.rho;
    rhoe = BL.rho(imax,:);
    
else
    rho  = ones(ny,nx);
    rhoe = ones(1,nx);
end

% Minus signs are due to the inverted mapping (freestream-to-wall)
BL.disp = -trapz(BL.y, 1-((rho.*BL.u)./repmat(rhoe.*BL.Ue,ny,1)) );
BL.mom  = -trapz(BL.y,   ((rho.*BL.u)./repmat(rhoe.*BL.Ue,ny,1)).*(1-(BL.u./repmat(BL.Ue,ny,1))) );
% Shape factor
BL.H12  =  BL.disp./BL.mom;

if BL.model
    % Stanton number
	dT = -FD1d2o_uneven(repmat(BL.y,1,nx),BL.T);
    if isempty(BL.Tw)
        Tw = zeros(1,nx);
    else
        Tw = BL.Tw;
    end
	BL.Stanton = BL.k(end,:).*dT(end,:)./(BL.rhoe.*BL.Ue.*BL.cpe.*(BL.Te*(1 + 0.84*(BL.gamma-1)/2*BL.Minf^2) - Tw));
end

% The 99% boundary layer thickness [m]
BL.edge = zeros(1,nx);
for i = 1:nx
    % Find last index where u is increasing
    id = find( diff(BL.u(:,i)./BL.u(imax,i)) >= 0 ,1,'last') + 1;
    % Only select monotonically decreasing part
    % (interp1 requires increasing u, hence the minus signs)
    BL.edge(i) = interp1(-BL.u(id:end,i)./BL.u(imax,i) ,BL.y(id:end), -0.99);
end

end