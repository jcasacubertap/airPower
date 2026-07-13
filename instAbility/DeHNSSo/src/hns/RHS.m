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
    R  = 0*StabRes.R;
    is_compr = isfield(Opt,'compressible') && strcmpi(Opt.compressible,'on');
    nvar     = 4 + double(is_compr);
    n_bc     = 3 + double(is_compr);

    % Station 1: inflow
    if ismember(j, Mode.RunJ)
        R(1:nvar*ny, 1) = Opt.AF * StabRes.phiIC(j,:);
    end

    % Stations 2..nx: inhomogeneous BCs (zero by default)
    if j <= size(Bound.bcw,3) && (any(Bound.bcw(:,:,j),'all') || any(Bound.bct(:,:,j),'all'))
        i_vec = 2:Grid.nx;
        base  = nvar*ny*(i_vec - 1);
        for k = 1:n_bc
            if strcmpi(Opt.bc_top{k}, 'IH_DR')
                R(base + (k-1)*ny + 1) = Bound.bct(k, i_vec, j);  % freestream
            end
            R(base + k*ny) = Bound.bcw(k, i_vec, j);              % wall
        end
    end
end

