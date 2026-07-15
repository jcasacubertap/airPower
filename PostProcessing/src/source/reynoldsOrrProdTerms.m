function sRO = reynoldsOrrProdTerms(sBF, sPert, inp)
% REYNOLDSORRPRODTERMS  Reynolds-Orr production of perturbation kinetic energy.
%
%   sRO = reynoldsOrrProdTerms(sBF, sPert, inp)
%
% Production term of the Reynolds-Orr equation for each spanwise-periodic
% perturbation mode, integrated in the spanwise (z) direction:
%
%   P = - lambda_z * sum_{i,j}  < u'_i u'_j >_z  dU_i/dx_j ,
%
%   i in {u,v,w},  j in {x,y}   (the base flow is spanwise-homogeneous, d/dz=0),
%   < u'_i u'_j >_z = Re( q_i conj(q_j) )     (see NORMALIZATION note below),
%   lambda_z = 2*pi/beta                       (from integrating 1 dz over a period).
%
% Full 3-D: the crossflow terms w', dW/dx, dW/dy are included. (The 2021
% original zeroed them for a 2-D TS case — see reynoldsOrrProdTermsFun_legacy2021.m.)
% Also returns the tangential/normal decomposition I1..I4 (I1+I2+I3+I4 = P),
% ported exactly from that script apart from the factor 2 (see NORMALIZATION).
%
% INPUTS (from importData, loadMode='loadFields'):
%   sBF   base flow on the stability grid, Ny x Nx, rows free-stream(1,:)->wall(end,:):
%         .x .y (meshgrid), .u .v .w, gradients .ux .uy .vx .vy .wx .wy
%   sPert perturbation modes (full, amplitude-scaled fields from importData):
%         .u .v .w  (Nmode x Ny x Nx) complex perturbation (NOT normalized)
%         .beta     (1 x Nmode)       spanwise wavenumber of each mode
%   inp.ro (all optional):
%         .modeIdx    modes to process (default: all with beta>0; beta=0 MFD excluded)
%         .xLim,.yLim integration window in base-flow units (default: full domain)
%
% OUTPUT sRO:
%   .beta    (1 x M)         spanwise wavenumber of each processed mode
%   .P       (Ny x Nx x M)   production field per mode
%   .intP    (1 x M)         volume integral of P over the window (signed)
%   .intAbsP (1 x M)         volume integral of |P| over the window
%   .peakAmp (1 x M)         max|u'| per mode (amplitude sanity-check)
%   .I1..I4  (Ny x Nx x M)   tangential/normal production decomposition:
%                            I1 normal-normal, I2 tang-normal, I3 normal-tang,
%                            I4 tang-tang.  I1+I2+I3+I4 = P.
%   .I2a..I2f, .I4a..I4f     per-base-flow-gradient sub-terms of I2 and I4
%   .x .y                    the grid (Ny x Nx)

    % ---- validate inputs ----
    need = {'x','y','u','v','w','ux','uy','vx','vy','wx','wy'};
    miss = need(~isfield(sBF, need));
    if ~isempty(miss)
        error('reynoldsOrr:sBF', 'sBF missing field(s): %s', strjoin(miss, ', '));
    end
    if ~all(isfield(sPert, {'u','v','w','beta'}))
        error('reynoldsOrr:sPert', ...
              'sPert needs .u/.v/.w (Nmode x Ny x Nx, full fields) and .beta (1 x Nmode).');
    end
    [Ny, Nx] = size(sBF.u);

    % ---- options (inp.ro overrides these defaults) ----
    ro = struct('modeIdx', [], 'xLim', [], 'yLim', []);
    if isfield(inp, 'ro') && ~isempty(inp.ro)
        f = fieldnames(inp.ro);
        for i = 1:numel(f)
            if isfield(ro, f{i}), ro.(f{i}) = inp.ro.(f{i}); end
        end
    end

    % ---- select spanwise-periodic modes (beta=0 MFD excluded: lambda_z -> Inf) ----
    beta   = sPert.beta(:).';
    bscale = max(1, max(abs(beta)));
    if isempty(ro.modeIdx)
        modeIdx = find(abs(beta) > 1e-8 * bscale);
    else
        modeIdx = ro.modeIdx(:).';
    end
    if isempty(modeIdx)
        error('reynoldsOrr:noModes', ...
              'No spanwise-periodic modes (beta>0). betavec = [%s].', num2str(beta));
    end
    M = numel(modeIdx);

    % ---- production field per mode ----
    sRO.x = sBF.x;  sRO.y = sBF.y;
    sRO.beta    = beta(modeIdx);
    sRO.P       = zeros(Ny, Nx, M);
    sRO.peakAmp = zeros(1, M);
    % tangential/normal decomposition fields (ported from the 2021 script)
    Idecomp = {'I1','I2','I3','I4', 'I2a','I2b','I2c','I2d','I2e','I2f', ...
               'I4a','I4b','I4c','I4d','I4e','I4f'};
    for i = 1:numel(Idecomp); sRO.(Idecomp{i}) = zeros(Ny,Nx,M); end

    for m = 1:M
        k   = modeIdx(m);
        lam = 2*pi / beta(k);                        % spanwise wavelength lambda_z

        % perturbation fields q_i(y,x) (full; DeHNSSo already combines mode + cc)
        qu = reshape(sPert.u(k,:,:), [Ny, Nx]);
        qv = reshape(sPert.v(k,:,:), [Ny, Nx]);
        qw = reshape(sPert.w(k,:,:), [Ny, Nx]);

        % z-averaged Reynolds stresses  < u'_i u'_j >_z = Re(q_i conj(q_j))  (no factor 2)
        Ruu = real(qu .* conj(qu));
        Ruv = real(qu .* conj(qv));                  % = Rvu
        Rvv = real(qv .* conj(qv));
        Rwu = real(qw .* conj(qu));                  % = Ruw
        Rwv = real(qw .* conj(qv));                  % = Rvw

        % production (z-integrated): P = -lambda_z * sum_ij <u'_i u'_j> dU_i/dx_j
        sRO.P(:,:,m) = -lam .* ( Ruu.*sBF.ux + Ruv.*sBF.uy ...
                               + Ruv.*sBF.vx + Rvv.*sBF.vy ...
                               + Rwu.*sBF.wx + Rwv.*sBF.wy );
        sRO.peakAmp(m) = max(abs(qu(:)));

        % ---- tangential/normal decomposition I1..I4 (ported from the 2021
        %      script; leading factor 2 dropped, as for P). Amplitude/phase form
        %      matching the original: A/2 = |q|, so (A/2)cos(Phase)=Re(q), etc.
        Au = 2*abs(qu);  Pu = angle(qu);
        Av = 2*abs(qv);  Pv = angle(qv);
        Aw = 2*abs(qw);  Pw = angle(qw);

        moduleBF = sqrt(sBF.u.^2 + sBF.v.^2 + sBF.w.^2);
        gamma1 = Au./2 .* sBF.u .* cos(Pu) + Av./2 .* sBF.v .* cos(Pv) + Aw./2 .* sBF.w .* cos(Pw);
        gamma2 = -Au./2 .* sBF.u .* sin(Pu) - Av./2 .* sBF.v .* sin(Pv) - Aw./2 .* sBF.w .* sin(Pw);
        gmag   = sqrt(gamma1.^2 + gamma2.^2);
        xi     = gmag ./ moduleBF.^2;

        % tangential vector
        shapeTangA = gmag ./ moduleBF.^2 .* sBF.u;
        shapeTangB = gmag ./ moduleBF.^2 .* sBF.v;
        shapeTangC = gmag ./ moduleBF.^2 .* sBF.w;
        phaseTang  = atan2(-gamma2, gamma1);
        FCTangA = shapeTangA .* exp(1i*phaseTang);
        FCTangB = shapeTangB .* exp(1i*phaseTang);
        FCTangC = shapeTangC .* exp(1i*phaseTang);

        % normal vector
        phaseNormA = atan2(Au./2.*sin(Pu) - sBF.u.*xi.*sin(phaseTang), Au./2.*cos(Pu) - sBF.u.*xi.*cos(phaseTang));
        phaseNormB = atan2(Av./2.*sin(Pv) - sBF.v.*xi.*sin(phaseTang), Av./2.*cos(Pv) - sBF.v.*xi.*cos(phaseTang));
        phaseNormC = atan2(Aw./2.*sin(Pw) - sBF.w.*xi.*sin(phaseTang), Aw./2.*cos(Pw) - sBF.w.*xi.*cos(phaseTang));
        shapeNormA = sqrt((Au./2).^2 + sBF.u.^2.*xi.^2 - 2.*sBF.u.*xi.*abs(Au./2).*cos(Pu - phaseTang));
        shapeNormB = sqrt((Av./2).^2 + sBF.v.^2.*xi.^2 - 2.*sBF.v.*xi.*abs(Av./2).*cos(Pv - phaseTang));
        shapeNormC = sqrt((Aw./2).^2 + sBF.w.^2.*xi.^2 - 2.*sBF.w.*xi.*abs(Aw./2).*cos(Pw - phaseTang));
        FCNormA = shapeNormA .* exp(1i*phaseNormA);
        FCNormB = shapeNormB .* exp(1i*phaseNormB);
        FCNormC = shapeNormC .* exp(1i*phaseNormC);

        % I1 (normal-normal)
        I1a = sBF.ux.*real(FCNormA.*conj(FCNormA));
        I1b = sBF.uy.*real(FCNormA.*conj(FCNormB));
        I1c = sBF.vx.*real(FCNormB.*conj(FCNormA));
        I1d = sBF.vy.*real(FCNormB.*conj(FCNormB));
        I1e = sBF.wx.*real(FCNormC.*conj(FCNormA));
        I1f = sBF.wy.*real(FCNormC.*conj(FCNormB));
        sRO.I1(:,:,m) = -lam .* (I1a+I1b+I1c+I1d+I1e+I1f);

        % I2 (tangential-normal)
        I2a = sBF.ux.*real(FCTangA.*conj(FCNormA));
        I2b = sBF.uy.*real(FCTangA.*conj(FCNormB));
        I2c = sBF.vx.*real(FCTangB.*conj(FCNormA));
        I2d = sBF.vy.*real(FCTangB.*conj(FCNormB));
        I2e = sBF.wx.*real(FCTangC.*conj(FCNormA));
        I2f = sBF.wy.*real(FCTangC.*conj(FCNormB));
        sRO.I2(:,:,m)  = -lam .* (I2a+I2b+I2c+I2d+I2e+I2f);
        sRO.I2a(:,:,m) = -lam.*I2a;  sRO.I2b(:,:,m) = -lam.*I2b;  sRO.I2c(:,:,m) = -lam.*I2c;
        sRO.I2d(:,:,m) = -lam.*I2d;  sRO.I2e(:,:,m) = -lam.*I2e;  sRO.I2f(:,:,m) = -lam.*I2f;

        % I3 (normal-tangential)
        I3a = sBF.ux.*real(FCNormA.*conj(FCTangA));
        I3b = sBF.uy.*real(FCNormA.*conj(FCTangB));
        I3c = sBF.vx.*real(FCNormB.*conj(FCTangA));
        I3d = sBF.vy.*real(FCNormB.*conj(FCTangB));
        I3e = sBF.wx.*real(FCNormC.*conj(FCTangA));
        I3f = sBF.wy.*real(FCNormC.*conj(FCTangB));
        sRO.I3(:,:,m) = -lam .* (I3a+I3b+I3c+I3d+I3e+I3f);

        % I4 (tangential-tangential)
        I4a = sBF.ux.*real(FCTangA.*conj(FCTangA));
        I4b = sBF.uy.*real(FCTangA.*conj(FCTangB));
        I4c = sBF.vx.*real(FCTangB.*conj(FCTangA));
        I4d = sBF.vy.*real(FCTangB.*conj(FCTangB));
        I4e = sBF.wx.*real(FCTangC.*conj(FCTangA));
        I4f = sBF.wy.*real(FCTangC.*conj(FCTangB));
        sRO.I4(:,:,m)  = -lam .* (I4a+I4b+I4c+I4d+I4e+I4f);
        sRO.I4a(:,:,m) = -lam.*I4a;  sRO.I4b(:,:,m) = -lam.*I4b;  sRO.I4c(:,:,m) = -lam.*I4c;
        sRO.I4d(:,:,m) = -lam.*I4d;  sRO.I4e(:,:,m) = -lam.*I4e;  sRO.I4f(:,:,m) = -lam.*I4f;
    end

    % Wall row (no-slip -> |V_base|=0 -> 0/0 in the tang/normal decomposition).
    % The 2021 script used cell-centre data with no wall row; here the wall is a
    % grid row, so zero it (the perturbation, hence its production, vanishes there).
    for i = 1:numel(Idecomp); sRO.(Idecomp{i})(end,:,:) = 0; end

    fprintf('reynoldsOrr: %d mode(s) [beta = %s]\n', M, num2str(sRO.beta, '%.4g '));

    % ---- volume integral of P over the (optional) window ----
    xrow  = sBF.x(1, :).';           % Nx x 1  (x varies along columns)
    ycol  = sBF.y(:, 1);             % Ny x 1  (y varies down rows)
    xmask = true(Nx,1);  if ~isempty(ro.xLim), xmask = xrow>=min(ro.xLim) & xrow<=max(ro.xLim); end
    ymask = true(Ny,1);  if ~isempty(ro.yLim), ymask = ycol>=min(ro.yLim) & ycol<=max(ro.yLim); end

    sRO.intP    = zeros(1, M);
    sRO.intAbsP = zeros(1, M);
    if nnz(xmask) < 2 || nnz(ymask) < 2
        warning('reynoldsOrr:window', ...
                'Integration window has <2 points in x or y; integrals left as 0.');
        return;
    end
    xi       = xrow(xmask);
    [yi, yo] = sort(ycol(ymask), 'ascend');          % ascending y so dy > 0 (physical)
    for m = 1:M
        Pm = sRO.P(ymask, xmask, m);
        Pm = Pm(yo, :);                              % rows -> ascending y
        sRO.intP(m)    = trapz(yi, trapz(xi, Pm,      2));
        sRO.intAbsP(m) = trapz(yi, trapz(xi, abs(Pm), 2));
    end
end
