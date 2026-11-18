function [Mode,StabRes] = inflow(Mode,Stab,StabGrid,StabRes,Opt)
%INFLOW  Set inflow boundary conditions for all active modes.
%
%   [Mode, StabRes] = inflow(Mode, Stab, StabGrid, StabRes)
%
%   For each mode with nonzero initial amplitude (Stab.A0), computes the
%   inflow condition using ILST (local stability eigenvalue problem) or
%   applies wall/loaded BCs. Sets phase reference (wall pressure or max u).
%   Populates StabRes.u/v/w/p/phi at station 1 and updates Mode.RunJ.

nonzeromodeswall = [];
nonzeromodesic   = [];

% Modes activated by wall/freestream actuation (Stab.bcw / Stab.bct)
nf_bc = min(size(Stab.bcw,3), Mode.nf);
for j = 1:nf_bc
    if any(Stab.bcw(:,:,j),'all') || any(Stab.bct(:,:,j),'all')
        nonzeromodeswall(end+1) = j; %#ok<AGROW>
    end
end

% Initialize solution vector phi = [u v w p]', Amplitude A and streamwise wavenumber alpha 
StabRes.phi    = zeros(Mode.nf,4*StabGrid.ny,StabGrid.nx);
StabRes.alpha  = zeros(Mode.nf,StabGrid.nx);
StabRes.A      = zeros(Mode.nf,StabGrid.nx);
StabRes.u      = zeros(Mode.nf,StabGrid.ny,StabGrid.nx);
StabRes.v      = zeros(Mode.nf,StabGrid.ny,StabGrid.nx);
StabRes.w      = zeros(Mode.nf,StabGrid.ny,StabGrid.nx);
StabRes.p      = zeros(Mode.nf,StabGrid.ny,StabGrid.nx);

% Determine ILST freestream BC from per-component Opt.bc_top.
% ILST uses a single flag: 'H_DR' or 'H_NM'.
% If any component is Neumann, use Neumann for ILST.
if ~exist('Opt','var') || ~isfield(Opt,'bc_top')
    bc_top = 'H_DR';
elseif iscell(Opt.bc_top)
    if any(strcmpi(Opt.bc_top, 'H_NM'))
        bc_top = 'H_NM';
    else
        bc_top = 'H_DR';
    end
else
    bc_top = Opt.bc_top;
end

if Stab.IC=="LOAD"
    % Load inflow condition from a previously saved StabRes or phiIC file.
    % The file must be specified via Stab.ICfile and contain either:
    %   - a variable 'StabRes' with fields u,v,w,p of size [nf x ny x nx], or
    %   - a variable 'phiIC'   of size [nf x 4*ny] (concatenated inflow eigenvectors).
    %
    % The loaded result must have matching ny and nf (same M, N, omega_0, beta_0).

    assert(isfield(Stab,'ICfile') && ~isempty(Stab.ICfile), ...
        'Stab.ICfile must be set when Stab.IC = ''LOAD''.');
    assert(isfile(Stab.ICfile), 'ICfile not found: %s', Stab.ICfile);

    loaded = load(Stab.ICfile);

    if isfield(loaded,'phiIC')
        % Direct load of phiIC [nf x 4*ny]
        phiIC_load = loaded.phiIC;
        assert(size(phiIC_load,1) == Mode.nf && size(phiIC_load,2) == 4*StabGrid.ny, ...
            'Loaded phiIC size [%d x %d] does not match expected [%d x %d].', ...
            size(phiIC_load,1), size(phiIC_load,2), Mode.nf, 4*StabGrid.ny);
        ny = StabGrid.ny;
        for j = 1:Mode.nf
            StabRes.u(j,:,1) = phiIC_load(j,      1:ny);
            StabRes.v(j,:,1) = phiIC_load(j,   ny+1:2*ny);
            StabRes.w(j,:,1) = phiIC_load(j, 2*ny+1:3*ny);
            StabRes.p(j,:,1) = phiIC_load(j, 3*ny+1:4*ny);

            %Hard code antonis
            yy = StabGrid.y;
            StabRes.u(j,:,1) = 0*phiIC_load(j,      1:ny);
            StabRes.v(j,:,1) = 1-yy.^2;
            StabRes.w(j,:,1) = abs(sin(pi*yy)) ;
            StabRes.p(j,:,1) = 0*phiIC_load(j, 3*ny+1:4*ny);




        end

    elseif isfield(loaded,'StabRes')
        SR = loaded.StabRes;
        assert(size(SR.u,1) == Mode.nf && size(SR.u,2) == StabGrid.ny, ...
            'Loaded StabRes.u size [%d x %d x ...] does not match expected [%d x %d x ...].', ...
            size(SR.u,1), size(SR.u,2), Mode.nf, StabGrid.ny);
        % Full-domain warm-start if the saved StabRes has the whole domain
        % (third dim matches StabGrid.nx). Otherwise fall back to station-1
        % only (original behaviour — inflow BC only).
        if size(SR.u,3) == StabGrid.nx
            % Amplitude scaling for continuation: if the saved run used a
            % different Stab.A0_fund, rescale the warm-start linearly so the
            % downstream saturation level is consistent with the new A0. The
            % nonlinear fixed point will shift slightly from this estimate, but
            % it's a far better initial guess than the raw saved values.
            scale = 1;
            if isfield(loaded,'Stab') && isfield(loaded.Stab,'A0_fund') ...
                    && isfield(Stab,'A0_fund') && loaded.Stab.A0_fund > 0
                scale = Stab.A0_fund / loaded.Stab.A0_fund;
            end
            for j = 1:Mode.nf
                StabRes.u(j,:,:) = scale * SR.u(j,:,:);
                StabRes.v(j,:,:) = scale * SR.v(j,:,:);
                StabRes.w(j,:,:) = scale * SR.w(j,:,:);
                StabRes.p(j,:,:) = scale * SR.p(j,:,:);
            end
            % Also seed StabRes.q (the flat solution vector Picard iterates on)
            % so the first Picard step sees the warm-start, not zeros.
            if isfield(SR,'q') && isequal(size(SR.q), [Mode.nf, 4*StabGrid.ny*StabGrid.nx])
                StabRes.q = scale * SR.q;
            end
            fprintf('  LOAD: full-domain warm-start (%d x %d x %d), A0 scale = %.3f\n', ...
                    size(SR.u,1), size(SR.u,2), size(SR.u,3), scale);
        else
            for j = 1:Mode.nf
                StabRes.u(j,:,1) = SR.u(j,:,1);
                StabRes.v(j,:,1) = SR.v(j,:,1);
                StabRes.w(j,:,1) = SR.w(j,:,1);
                StabRes.p(j,:,1) = SR.p(j,:,1);
            end
            fprintf('  LOAD: inflow-only (station 1) — no downstream warm-start\n');
        end

    else
        error('ICfile must contain a variable named ''StabRes'' or ''phiIC''.');
    end

    % Apply amplitude scaling and phase reference for physical modes
    [~, nonzeromodesic] = find(Stab.A0 > 0);
    for j = flip(nonzeromodesic)
        if j >= round(Mode.nf/2)
            % Normalise by max(|u|) before scaling by A0
            y_int  = linspace(StabGrid.y(1), StabGrid.y(end), 4000);
            u1int  = interp1(StabGrid.y, abs(StabRes.u(j,:,1)), y_int, 'spline');
            ampltd = max(abs(u1int));
            if ampltd > 0
                StabRes.u(j,:,1) = StabRes.u(j,:,1) / ampltd;
                StabRes.v(j,:,1) = StabRes.v(j,:,1) / ampltd;
                StabRes.w(j,:,1) = StabRes.w(j,:,1) / ampltd;
                StabRes.p(j,:,1) = StabRes.p(j,:,1) / ampltd;
            end
            StabRes.u(j,:,1) = Stab.A0(j) .* StabRes.u(j,:,1);
            StabRes.v(j,:,1) = Stab.A0(j) .* StabRes.v(j,:,1);
            StabRes.w(j,:,1) = Stab.A0(j) .* StabRes.w(j,:,1);
            StabRes.p(j,:,1) = Stab.A0(j) .* StabRes.p(j,:,1);
            StabRes.phi(j,:,1) = [StabRes.u(j,:,1) StabRes.v(j,:,1) StabRes.w(j,:,1) StabRes.p(j,:,1)];
            switch lower(Stab.phaseRef)
                case 'pwall'
                    theta0 = angle(StabRes.p(j,end,1));
                case 'umax'
                    [~,iy] = max(abs(StabRes.u(j,:,1)));
                    theta0 = angle(StabRes.u(j,iy,1));
                otherwise
                    error('Unknown phaseRef option');
            end
            StabRes.phi(j,:,1) = StabRes.phi(j,:,1) * exp(-1i*theta0);
        else
            % Conjugate mode
            j2 = round(Mode.nf/2) - (j - round(Mode.nf/2));
            StabRes.u(j,:,1)   = conj(StabRes.u(j2,:,1));
            StabRes.v(j,:,1)   = conj(StabRes.v(j2,:,1));
            StabRes.w(j,:,1)   = conj(StabRes.w(j2,:,1));
            StabRes.p(j,:,1)   = conj(StabRes.p(j2,:,1));
            StabRes.alpha(j,1) = -conj(StabRes.alpha(j2,1));
            StabRes.phi(j,:,1) = conj(StabRes.phi(j2,:,1));
        end
    end
    fprintf('Initialising with LOAD solution from: %s\n', Stab.ICfile);

elseif Stab.IC=="ILST"
    [~,nonzeromodesic]= find(Stab.A0>0);

    for j = flip(nonzeromodesic) %Run initial condition for all nonzero modes
    if j >= round(Mode.nf/2)% Mode is physical
        [StabRes]=ILST(j,1,StabGrid,StabRes,bc_top,Opt);
    
        % Output from IC_HNS is normalized by max(u'), superimpose amplitude
        StabRes.u(j,:,1) = Stab.A0(j).*StabRes.u(j,:,1);
        StabRes.v(j,:,1) = Stab.A0(j).*StabRes.v(j,:,1);
        StabRes.w(j,:,1) = Stab.A0(j).*StabRes.w(j,:,1);
        StabRes.p(j,:,1) = Stab.A0(j).*StabRes.p(j,:,1);
        StabRes.phi(j,:,1) = [StabRes.u(j,:,1) StabRes.v(j,:,1) StabRes.w(j,:,1) StabRes.p(j,:,1)];
        switch lower(Stab.phaseRef)  % e.g. 'pwall' or 'umax'
        case 'pwall'
            theta0 = angle(StabRes.p(j,end,1));  % wall pressure
        case 'umax'
            [~,iy] = max(abs(StabRes.u(j,:,1))); % max |u'|
            theta0 = angle(StabRes.u(j,iy,1));
        otherwise
            error('Unknown phaseRef option');
        end
        % reference the phase of the inflow perturbation
        StabRes.phi(j,:,1) = StabRes.phi(j,:,1) * exp(-1i*theta0);
    
    else % Apply symmetry for conjugate mode
        j2 = round(Mode.nf/2)-(j-round(Mode.nf/2));% Defines conjugate mode
        
        StabRes.u(j,:,1) = conj(StabRes.u(j2,:,1));
        StabRes.v(j,:,1) = conj(StabRes.v(j2,:,1));
        StabRes.w(j,:,1) = conj(StabRes.w(j2,:,1));
        StabRes.p(j,:,1) = conj(StabRes.p(j2,:,1));
        StabRes.alpha(j,1) = -conj(StabRes.alpha(j,1)); % growth rates should be equal between conjugate modes so -conj
        StabRes.phi(j,:,1) = conj(StabRes.phi(j2,:,1));
    end                                                                            
    end 
    
    if max(abs(Stab.bcw(2:end,:)))==0
        fprintf('Initialising with ILST solution at inflow. \n')
    else
        fprintf('Initialising with ILST solution at inflow and non-homogeneous b.c. at wall. \n')
    end

else
    error('Unknown inflow condition Stab.IC = ''%s''. Valid options: ''ILST'', ''LOAD''.', Stab.IC);
end

% nonzero modes are those initiated at inflow + those forced at wall
nonzeromodes = unique([nonzeromodesic, nonzeromodeswall]);
Mode.RunJ = nonzeromodes;            % Used in loops to assess only active modes
Mode.Aactive(Mode.RunJ) = 1;      % Set active modes to introduced modes

% Define phiIC
StabRes.phiIC = squeeze(StabRes.phi(:,:,1));

end

