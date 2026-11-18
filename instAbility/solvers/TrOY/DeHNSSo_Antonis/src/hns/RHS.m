function [R] = RHS(Bound,Grid,Mode,Opt,StabRes,j)
%RHS  Construct the right-hand side vector for mode j.
%
%   R = RHS(Bound, Grid, Mode, Opt, StabRes, j)
%
%   Assembles the RHS from:
%     - Station 1: inflow shape function (phiIC), scaled by amplitude factor AF.
%     - Stations 2..nx: inhomogeneous wall/freestream BCs from Bound.bcw/bct.
%       Grid convention: index 1 = freestream (y=H), index ny = wall (y=0).
%       Freestream values injected only for IH_DR components (Opt.bc_top).
%
%   Output: R [N x 1] — RHS vector for the global linear system ML * q = R.

    ny = Grid.ny;
    R = 0*StabRes.R;

    % Station 1: inflow
    if ismember(j, Mode.RunJ)
        R(1:4*ny, 1) = Opt.AF * StabRes.phiIC(j,:);
    end

    % Stations 2..nx: inhomogeneous BCs (zero by default, non-zero for actuation)
    if j <= size(Bound.bcw,3) && (any(Bound.bcw(:,:,j),'all') || any(Bound.bct(:,:,j),'all'))
        i_vec = 2:Grid.nx;
        base  = 4*ny*(i_vec - 1);   % 1 x (nx-1) base offsets
        for k = 1:3
            % Freestream: only inject for IH_DR (inhomogeneous Dirichlet)
            if strcmpi(Opt.bc_top{k}, 'IH_DR')
                R(base + (k-1)*ny + 1) = Bound.bct(k, i_vec, j);  % freestream
            end
            % Wall: always inject (wall BCs are always Dirichlet)
            R(base + k*ny) = Bound.bcw(k, i_vec, j);              % wall
        end
    end
end

