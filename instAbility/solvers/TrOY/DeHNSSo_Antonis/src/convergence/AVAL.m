function [Opt,StabRes] = AVAL(Grid,Opt,Mode,StabRes)
%AVAL  Amplitude value check and inflow amplitude damping.
%
%   [Opt, StabRes] = AVAL(Grid, Opt, Mode, StabRes)
%
%   Checks if the maximum perturbation amplitude exceeds Opt.AMAX.
%   If so, reduces the inflow amplitude by half and resets the solution.
%   Tracks the number of damping events in Opt.AVAL.

    % Reshape q into phi (vectorised)
    ny_ = Grid.ny;
    nx_ = Grid.nx;
    for j = Mode.RunJ
        StabRes.phi(j,:,:) = reshape(StabRes.q(j,:), 4*ny_, nx_);
    end

    % Calculate the maximum amplitude reached linearly
    [Opt.AVAL,~] = max(max(abs(squeeze(StabRes.phi(round(Mode.nf/2)+1,1:3*ny_,:)))));

    % Find Amplitude Factor (AF) such that evolution is linear
    if Opt.AVAL >= Opt.AMAX
        Opt.AF = Opt.AMAX/Opt.AVAL;
    end

    StabRes.q = Opt.AF*StabRes.q;
end
