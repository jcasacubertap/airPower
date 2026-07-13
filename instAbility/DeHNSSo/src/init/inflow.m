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

is_compr = isfield(Opt,'compressible') && strcmpi(Opt.compressible,'on');
nvar     = 4 + double(is_compr);

StabRes.phi    = zeros(Mode.nf,nvar*StabGrid.ny,StabGrid.nx);
StabRes.alpha  = zeros(Mode.nf,StabGrid.nx);
StabRes.A      = zeros(Mode.nf,StabGrid.nx);
StabRes.u      = zeros(Mode.nf,StabGrid.ny,StabGrid.nx);
StabRes.v      = zeros(Mode.nf,StabGrid.ny,StabGrid.nx);
StabRes.w      = zeros(Mode.nf,StabGrid.ny,StabGrid.nx);
StabRes.p      = zeros(Mode.nf,StabGrid.ny,StabGrid.nx);
if is_compr
    StabRes.rho = zeros(Mode.nf,StabGrid.ny,StabGrid.nx);
    StabRes.T   = zeros(Mode.nf,StabGrid.ny,StabGrid.nx);
end

if is_compr && ~(strcmpi(Stab.IC,'CLST') || strcmpi(Stab.IC,'LOAD'))
    error('inflow:badIC', ...
        'Opt.compressible=''on'' requires Stab.IC=''CLST'' or ''LOAD''; got ''%s''.', Stab.IC);
end
if ~is_compr && strcmpi(Stab.IC,'CLST')
    error('inflow:badIC', 'Stab.IC=''CLST'' requires Opt.compressible=''on''.');
end

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
    % phiIC width = nvar*ny: 4 incompressible (u,v,w,p), 5 compressible (u,v,w,ρ,T).

    assert(isfield(Stab,'ICfile') && ~isempty(Stab.ICfile), ...
        'Stab.ICfile must be set when Stab.IC = ''LOAD''.');
    assert(isfile(Stab.ICfile), 'ICfile not found: %s', Stab.ICfile);

    loaded = load(Stab.ICfile);

    if isfield(loaded,'phiIC')
        phiIC_load = loaded.phiIC;
        assert(size(phiIC_load,1) == Mode.nf && size(phiIC_load,2) == nvar*StabGrid.ny, ...
            'Loaded phiIC size [%d x %d] does not match expected [%d x %d].', ...
            size(phiIC_load,1), size(phiIC_load,2), Mode.nf, nvar*StabGrid.ny);
        ny = StabGrid.ny;
        for j = 1:Mode.nf
            StabRes.u(j,:,1) = phiIC_load(j,      1:ny);
            StabRes.v(j,:,1) = phiIC_load(j,   ny+1:2*ny);
            StabRes.w(j,:,1) = phiIC_load(j, 2*ny+1:3*ny);
            if is_compr
                StabRes.rho(j,:,1) = phiIC_load(j, 3*ny+1:4*ny);
                StabRes.T(j,:,1)   = phiIC_load(j, 4*ny+1:5*ny);
            else
                StabRes.p(j,:,1)   = phiIC_load(j, 3*ny+1:4*ny);
            end
        end

    elseif isfield(loaded,'StabRes')
        SR = loaded.StabRes;
        assert(size(SR.u,1) == Mode.nf && size(SR.u,2) == StabGrid.ny, ...
            'Loaded StabRes.u size [%d x %d x ...] does not match expected [%d x %d x ...].', ...
            size(SR.u,1), size(SR.u,2), Mode.nf, StabGrid.ny);
        if is_compr
            assert(isfield(SR,'rho') && isfield(SR,'T'), ...
                'Loaded StabRes lacks .rho/.T fields (compressible LOAD requires them).');
        end
        % Continuation: rescale the warm-start linearly when Stab.A0_fund changed.
        scale = 1;
        if isfield(loaded,'Stab') && isfield(loaded.Stab,'A0_fund') ...
                && isfield(Stab,'A0_fund') && loaded.Stab.A0_fund > 0
            scale = Stab.A0_fund / loaded.Stab.A0_fund;
        end
        if size(SR.u,3) == StabGrid.nx
            for j = 1:Mode.nf
                StabRes.u(j,:,:) = scale * SR.u(j,:,:);
                StabRes.v(j,:,:) = scale * SR.v(j,:,:);
                StabRes.w(j,:,:) = scale * SR.w(j,:,:);
                if is_compr
                    StabRes.rho(j,:,:) = scale * SR.rho(j,:,:);
                    StabRes.T(j,:,:)   = scale * SR.T(j,:,:);
                else
                    StabRes.p(j,:,:) = scale * SR.p(j,:,:);
                end
            end
            if isfield(SR,'q') && isequal(size(SR.q), [Mode.nf, nvar*StabGrid.ny*StabGrid.nx])
                StabRes.q = scale * SR.q;
            end
            fprintf('  LOAD: full-domain warm-start (%d x %d x %d), A0 scale = %.3f\n', ...
                    size(SR.u,1), size(SR.u,2), size(SR.u,3), scale);
        else
            for j = 1:Mode.nf
                StabRes.u(j,:,1) = SR.u(j,:,1);
                StabRes.v(j,:,1) = SR.v(j,:,1);
                StabRes.w(j,:,1) = SR.w(j,:,1);
                if is_compr
                    StabRes.rho(j,:,1) = SR.rho(j,:,1);
                    StabRes.T(j,:,1)   = SR.T(j,:,1);
                else
                    StabRes.p(j,:,1) = SR.p(j,:,1);
                end
            end
            fprintf('  LOAD: inflow-only (station 1) — no downstream warm-start\n');
        end

    else
        error('ICfile must contain a variable named ''StabRes'' or ''phiIC''.');
    end

    if is_compr
        % Derive p' from EOS:  p' = (ρ'·T̄ + ρ̄·T') / (γ·Ma²)
        gM2 = StabGrid.gamma * StabGrid.Ma^2;
        for j = 1:Mode.nf
            StabRes.p(j,:,1) = (StabRes.rho(j,:,1) .* StabGrid.T(:,1).' ...
                              + StabGrid.rho(:,1).' .* StabRes.T(j,:,1)) / gM2;
        end
    end

    [~, nonzeromodesic] = find(Stab.A0 > 0);
    for j = flip(nonzeromodesic)
        if j >= round(Mode.nf/2)
            % Normalise by max(|u|), then scale by A0.
            y_int  = linspace(StabGrid.y(1), StabGrid.y(end), 4000);
            u1int  = interp1(StabGrid.y, abs(StabRes.u(j,:,1)), y_int, 'spline');
            ampltd = max(abs(u1int));
            if ampltd > 0
                StabRes.u(j,:,1) = StabRes.u(j,:,1) / ampltd;
                StabRes.v(j,:,1) = StabRes.v(j,:,1) / ampltd;
                StabRes.w(j,:,1) = StabRes.w(j,:,1) / ampltd;
                if is_compr
                    StabRes.rho(j,:,1) = StabRes.rho(j,:,1) / ampltd;
                    StabRes.T(j,:,1)   = StabRes.T(j,:,1)   / ampltd;
                else
                    StabRes.p(j,:,1) = StabRes.p(j,:,1) / ampltd;
                end
            end
            StabRes.u(j,:,1) = Stab.A0(j) .* StabRes.u(j,:,1);
            StabRes.v(j,:,1) = Stab.A0(j) .* StabRes.v(j,:,1);
            StabRes.w(j,:,1) = Stab.A0(j) .* StabRes.w(j,:,1);
            if is_compr
                StabRes.rho(j,:,1) = Stab.A0(j) .* StabRes.rho(j,:,1);
                StabRes.T(j,:,1)   = Stab.A0(j) .* StabRes.T(j,:,1);
                StabRes.phi(j,:,1) = [StabRes.u(j,:,1) StabRes.v(j,:,1) StabRes.w(j,:,1) StabRes.rho(j,:,1) StabRes.T(j,:,1)];
            else
                StabRes.p(j,:,1) = Stab.A0(j) .* StabRes.p(j,:,1);
                StabRes.phi(j,:,1) = [StabRes.u(j,:,1) StabRes.v(j,:,1) StabRes.w(j,:,1) StabRes.p(j,:,1)];
            end
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
            if is_compr
                StabRes.rho(j,:,1) = conj(StabRes.rho(j2,:,1));
                StabRes.T(j,:,1)   = conj(StabRes.T(j2,:,1));
            else
                StabRes.p(j,:,1)   = conj(StabRes.p(j2,:,1));
            end
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

elseif Stab.IC=="CLST"
    [~,nonzeromodesic] = find(Stab.A0>0);

    for j = flip(nonzeromodesic)
        if j >= round(Mode.nf/2)   % Physical mode
            StabRes = CLST(j, 1, StabGrid, StabRes, bc_top, Opt);

            % CLST normalizes by max|u'|; scale by A0 and rebuild phi at the new amplitude.
            StabRes.u(j,:,1)   = Stab.A0(j).*StabRes.u(j,:,1);
            StabRes.v(j,:,1)   = Stab.A0(j).*StabRes.v(j,:,1);
            StabRes.w(j,:,1)   = Stab.A0(j).*StabRes.w(j,:,1);
            StabRes.p(j,:,1)   = Stab.A0(j).*StabRes.p(j,:,1);
            StabRes.rho(j,:,1) = Stab.A0(j).*StabRes.rho(j,:,1);
            StabRes.T(j,:,1)   = Stab.A0(j).*StabRes.T(j,:,1);

            StabRes.phi(j,:,1) = [StabRes.u(j,:,1),    StabRes.v(j,:,1), ...
                                  StabRes.w(j,:,1),    StabRes.rho(j,:,1), ...
                                  StabRes.T(j,:,1)];

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
        else   % Conjugate mode
            j2 = round(Mode.nf/2)-(j-round(Mode.nf/2));
            StabRes.u(j,:,1)   = conj(StabRes.u(j2,:,1));
            StabRes.v(j,:,1)   = conj(StabRes.v(j2,:,1));
            StabRes.w(j,:,1)   = conj(StabRes.w(j2,:,1));
            StabRes.p(j,:,1)   = conj(StabRes.p(j2,:,1));
            StabRes.rho(j,:,1) = conj(StabRes.rho(j2,:,1));
            StabRes.T(j,:,1)   = conj(StabRes.T(j2,:,1));
            StabRes.alpha(j,1) = -conj(StabRes.alpha(j2,1));
            StabRes.phi(j,:,1) = conj(StabRes.phi(j2,:,1));
        end
    end
    fprintf('Initialising with CLST solution at inflow.\n')

else
    error('Unknown inflow condition Stab.IC = ''%s''. Valid options: ''ILST'', ''CLST'', ''LOAD''.', Stab.IC);
end

% nonzero modes are those initiated at inflow + those forced at wall
nonzeromodes = unique([nonzeromodesic, nonzeromodeswall]);
Mode.RunJ = nonzeromodes;            % Used in loops to assess only active modes
Mode.Aactive(Mode.RunJ) = 1;      % Set active modes to introduced modes

% Define phiIC
StabRes.phiIC = squeeze(StabRes.phi(:,:,1));

end

