function [sBF, sPert, inp] = importData(inp)
% importData  Load base flow (and optionally perturbation) for post-processing.
%
%   [sBF, sPert, inp] = importData(inp)
%
% Dispatches on inp.loadMode:
%   'loadBF'     -> base flow only, from a PreProcessing case's midPlane data.
%                   sPert is returned empty ([]).
%   'loadFields' -> base flow AND perturbation, from io/input/<inp.fieldsFile>.
%                   The base flow has been interpolated onto a stability grid,
%                   so it differs from the midPlane base flow of 'loadBF'.
%
% Returns sBF with fields (Ny x Nx matrices):
%   sBF.x, sBF.y           meshgrid coordinates
%   sBF.u, sBF.v, sBF.w    velocity components
%   sBF.p                  pressure          (loadBF; absent in loadFields)
%   sBF.omz                z-vorticity       (loadBF, if present in source)
%   sBF.ux,uy,vx,vy,wx,wy  base-flow gradients (loadFields: read from source;
%                          loadBF: computed here from u,v,w via differentiateField)
%
% sPert (loadFields only) carries all perturbation modes:
%   sPert.u,v,w,p          (Nmode x Ny x Nx) complex, PHYSICAL PEAK fields
%                          (A .* shape, A = max|u'|) — use directly; RMS = |.|/sqrt(2).
%   sPert.A                (Nmode x Nx) raw StabRes.A (= max|u'| peak), diagnostic
%   sPert.omega, sPert.beta, sPert.alpha   (if present in source)
%
% Required inp fields:
%   inp.airPowerRoot   absolute path to airPower
%   inp.loadMode       'loadBF' (default) or 'loadFields'
%   inp.caseType       'DFP' or 'TTCP'        (loadBF)
%   inp.fieldsFile     name of .mat in io/input (loadFields)

    if ~isfield(inp, 'loadMode') || isempty(inp.loadMode)
        inp.loadMode = 'loadBF';
    end

    sPert = [];
    switch inp.loadMode
        case 'loadBF'
            [sBF, inp] = importFromPreProc(inp);
        case 'loadFields'
            [sBF, sPert, inp] = importFromFields(inp);
        otherwise
            error('importData:badLoadMode', ...
                  'Unknown loadMode: ''%s''. Must be ''loadBF'' or ''loadFields''.', ...
                  inp.loadMode);
    end

    % Enforce the post-processing row convention on every assembled matrix:
    %   row 1   (1,:)   = free-stream
    %   row end (end,:) = wall
    % Applied here, at the single load entry point, so all downstream consumers
    % and any matrix assembly (e.g. reynoldsOrrProdTerms) can rely on it. The
    % base-flow gradients (ux..wy) are Ny x Nx too, so they flip together with
    % the velocity fields and stay co-registered.
    [sBF, sPert] = orientFreestreamToWall(sBF, sPert);
end

% --- enforce row convention: (1,:) free-stream, (end,:) wall ---
% The wall-normal orientation is decided physically: the free-stream carries the
% largest velocity magnitude and the wall the smallest (no-slip). This is
% independent of how each source happened to store the wall-normal coordinate
% (loadBF sorts y ascending; loadFields takes the .mat's order as-is), so both
% paths end up consistent. The same flip is applied to the base flow and to
% every perturbation mode so they remain co-registered on the grid.
function [sBF, sPert] = orientFreestreamToWall(sBF, sPert)

    if ~isfield(sBF, 'u') || isempty(sBF.u)
        return;   % no velocity to decide on; leave as-is
    end
    [Ny, Nx] = size(sBF.u);

    % per-row mean speed (averaged over x), from whatever components exist
    spd2 = zeros(Ny, Nx);
    for c = {'u', 'v', 'w'}
        if isfield(sBF, c{1}) && isequal(size(sBF.(c{1})), [Ny, Nx])
            spd2 = spd2 + sBF.(c{1}).^2;
        end
    end
    rowSpeed = mean(sqrt(spd2), 2, 'omitnan');   % Ny x 1

    % Need both extremes to decide; bail (no flip) if either is undefined.
    if isnan(rowSpeed(1)) || isnan(rowSpeed(end))
        warning('importData:cannotOrient', ...
                ['Could not determine free-stream/wall orientation (NaN in the ', ...
                 'edge rows). Leaving row order unchanged — verify (1,:) is the ', ...
                 'free-stream and (end,:) the wall.']);
        return;
    end

    if rowSpeed(1) >= rowSpeed(end)
        fprintf('importData: rows already free-stream (1,:) -> wall (end,:)\n');
        return;
    end

    % row 1 is the low-speed (wall) side -> flip the wall-normal order
    flds = fieldnames(sBF);
    for i = 1:numel(flds)
        if isequal(size(sBF.(flds{i})), [Ny, Nx])
            sBF.(flds{i}) = flipud(sBF.(flds{i}));
        end
    end
    if ~isempty(sPert)
        pf = fieldnames(sPert);
        for i = 1:numel(pf)
            v = sPert.(pf{i});
            if ndims(v) == 3 && size(v, 2) == Ny     % (Nmode x Ny x Nx) modes
                sPert.(pf{i}) = flip(v, 2);
            end
        end
    end
    fprintf('importData: flipped wall-normal rows -> free-stream (1,:), wall (end,:)\n');
end

% --- loadBF: read a PreProcessing case's midPlane and arrange it on a meshgrid ---
% Reads <case>/postProcessing/midPlane.{bin,csv} from the case selected by
% inp.caseType. If both .bin and .csv are present, .bin wins. Also sets
% inp.wallExtrap depending on whether the file contained wall-extrapolation
% rows (u=v=w=0 added at wall face centres by the writeMidPlane function object).
function [sBF, inp] = importFromPreProc(inp)

    caseDir = caseDirFor(inp);

    % --- read source: prefer binary, fall back to CSV ---
    csvPath = fullfile(caseDir, 'postProcessing', 'midPlane.csv');
    binPath = fullfile(caseDir, 'postProcessing', 'midPlane.bin');
    if isfile(binPath)
        T = readMidPlaneBinary(binPath);
        fprintf('importData: read %s (%d rows)\n', binPath, height(T));
    elseif isfile(csvPath)
        T = readtable(csvPath);
        fprintf('importData: read %s (%d rows)\n', csvPath, height(T));
    else
        error('importData:notFound', ...
              'midPlane.bin/.csv not found in %s/postProcessing/', caseDir);
    end

    % Absolute tolerance for grid-position dedup. The ASCII midPlane.csv
    % is written with controlDict's writePrecision (typically 6 sig figs),
    % so cell positions of O(0.1 m) carry ~1e-7 m of FP noise after
    % round-trip. 1e-7 m collapses that noise without losing real cells.
    tol = 1e-7;

    % --- detect whether the producer included wall-extrapolation rows ---
    % The writeMidPlane function object writes wall-face centres with
    % u=v=w=0 (literal 0.0) when inputs.jl has wallExtrapolation=true.
    % Reflect that fact back into inp for downstream consumers.
    isWall = (T.u == 0) & (T.v == 0) & (T.w == 0);
    inp.wallExtrap = any(isWall);
    fprintf('importData: wallExtrap = %s (%d wall rows of %d total)\n', ...
            mat2str(inp.wallExtrap), nnz(isWall), height(T));

    % --- safety: the writeMidPlane function object must emit exactly one
    %     z-plane. If we see more than one distinct z value in the data,
    %     something on the OpenFOAM side has regressed (e.g. the z-tolerance
    %     widened, or the mesh got multiple cells per plane). Fail, so the
    %     user fixes the producer (system/controlDict, writeMidPlane).
    if ismember('z', T.Properties.VariableNames)
        zu = uniquetol(T.z, tol, 'DataScale', 1);
        if numel(zu) > 1
            error('importData:multipleZPlanes', ...
                ['midPlane data contains %d distinct z values (%s).\n', ...
                 'Expected exactly one — check writeMidPlane in:\n  ', ...
                 '<case>/system/controlDict\n', ...
                 'The cell-collection block should pick a single z-plane ', ...
                 '(nearest cell-center to zTarget), not a tolerance window.'], ...
                numel(zu), strjoin(compose('%.6g', zu), ', '));
        end
    end

    % --- arrange on a 2D meshgrid ---
    [xu, ~, ix] = uniquetol(T.x, tol, 'DataScale', 1);
    [yu, ~, iy] = uniquetol(T.y, tol, 'DataScale', 1);
    Nx = numel(xu);
    Ny = numel(yu);

    if Nx * Ny ~= height(T)
        warning('importData:nonRectangular', ...
                ['midPlane has %d rows but unique x * unique y = %d * %d = %d. ', ...
                 'Grid is not fully populated; empty cells will be NaN.'], ...
                height(T), Nx, Ny, Nx * Ny);
    end

    [Xg, Yg] = meshgrid(xu, yu);   % Ny x Nx
    sBF.x = Xg;
    sBF.y = Yg;

    lin = sub2ind([Ny, Nx], iy, ix);

    fields = {'u', 'v', 'w', 'p', 'omz'};
    for k = 1:numel(fields)
        f = fields{k};
        if ~ismember(f, T.Properties.VariableNames), continue; end
        M = nan(Ny, Nx);
        M(lin) = T.(f);
        sBF.(f) = M;
    end

    % --- base-flow gradients on the imported grid -----------------------
    % Computed here so the loadBF base flow carries the same ux..wy fields as
    % the loadFields path. differentiateField is 2nd-order on the non-uniform
    % grid and coordinate-based (it differences sBF.x/sBF.y, not row/column
    % indices), so the later free-stream/wall row flip in importData reorders
    % these derivatives consistently with the velocity fields they come from.
    %   <c>x = d<c>/dx (across columns, dim 2);  <c>y = d<c>/dy (down rows, dim 1)
    if Nx >= 3 && Ny >= 3
        for c = {'u', 'v', 'w'}
            comp = c{1};
            if ~isfield(sBF, comp), continue; end
            sBF.([comp 'x']) = differentiateField(sBF.(comp), sBF.x, 2);
            sBF.([comp 'y']) = differentiateField(sBF.(comp), sBF.y, 1);
        end
    else
        warning('importData:gradientsSkipped', ...
                ['Grid too small for 2nd-order base-flow gradients ', ...
                 '(Ny=%d, Nx=%d); ux..wy not computed.'], Ny, Nx);
    end
end

% --- loadFields: read base flow + perturbation from an io/input/*.mat ---
% The file is expected to hold (Casacuberta2022.mat-style):
%   StabGrid : base flow on the stability grid
%              .x (1 x Nx), .y (1 x Ny), .U/.V/.W and gradients .dxU/.dyU/...
%              (Ny x Nx)
%   StabRes  : perturbation modes, .u/.v/.w/.p (Nmode x Ny x Nx) complex
function [sBF, sPert, inp] = importFromFields(inp)

    if ~isfield(inp, 'fieldsFile') || isempty(inp.fieldsFile)
        error('importData:noFieldsFile', ...
              'inp.fieldsFile is required for loadMode=''loadFields''.');
    end

    matPath = fullfile(inp.airPowerRoot, 'PostProcessing', 'io', 'input', inp.fieldsFile);
    if ~isfile(matPath)
        error('importData:fieldsNotFound', 'Fields file not found: %s', matPath);
    end

    S = load(matPath);
    if ~isfield(S, 'StabGrid')
        error('importData:noStabGrid', ...
              '%s has no StabGrid struct (base flow on the stability grid).', matPath);
    end
    if ~isfield(S, 'StabRes')
        error('importData:noStabRes', ...
              '%s has no StabRes struct (perturbation modes).', matPath);
    end
    G = S.StabGrid;
    R = S.StabRes;

    % --- base flow on the stability grid (Ny x Nx) ---
    % Coordinates: prefer the PHYSICAL Cartesian arrays G.X/G.Y. They carry the
    % actual wall geometry (e.g. a hump/roughness: G.Y rises with the bump) and
    % use the inlet-origin frame, so the base flow matches the PreProcessing view.
    % The 1-D G.x/G.y are wall-fitted (x from the LE; eta = flat wall at 0), which
    % would flatten any wall feature. Fall back to meshgrid(G.x,G.y) for older
    % grids without X/Y (flat wall, where the two representations coincide).
    [Ny, Nx] = size(G.U);
    if isfield(G,'X') && isfield(G,'Y') && ...
       isequal(size(G.X),[Ny,Nx]) && isequal(size(G.Y),[Ny,Nx])
        Xg = G.X;  Yg = G.Y;
    else
        [Xg, Yg] = meshgrid(G.x(:).', G.y(:).');   % Ny x Nx (flat-wall fallback)
    end
    if ~isequal(size(Xg), [Ny, Nx])
        error('importData:fieldsGridMismatch', ...
              ['StabGrid base flow is %dx%d but the coordinate grid is %dx%d. ', ...
               'Expected y along rows and x along cols.'], ...
              Ny, Nx, size(Xg,1), size(Xg,2));
    end

    sBF.x = Xg;
    sBF.y = Yg;
    sBF.u = G.U;
    sBF.v = G.V;
    sBF.w = G.W;
    % base-flow gradients (available because already on the stability grid)
    sBF.ux = G.dxU;  sBF.uy = G.dyU;
    sBF.vx = G.dxV;  sBF.vy = G.dyV;
    sBF.wx = G.dxW;  sBF.wy = G.dyW;

    % --- perturbation: all modes (mode selection deferred downstream) ---
    % Store the PHYSICAL PEAK fields so u/v/w/p ARE the perturbation and every
    % downstream script uses them directly (no shape/A split). StabRes.u/v/w/p are
    % peak-normalized shapes and R.A the modal amplitude (Nmode x Nx). Per DeHNSSo
    % (postprocess.m / plot_amplitudes.m), StabRes.A = 2|u_max| = max|u'|, i.e. the
    % PEAK of the physical real perturbation  u'(z) = 2|a|cos(bz) = A*shape*cos.
    % So we fold in the FULL A (the peak) and the spanwise RMS is |sPert|/sqrt(2)
    % (DeHNSSo's own 'urms' metric = A/sqrt(2)). The complex Fourier coefficient is
    % A/2*shape; do NOT confuse it with the peak.
    if isfield(R, 'A')
        Am = reshape(R.A, size(R.A,1), 1, size(R.A,2));   % Nmode x 1 x Nx (full A = peak, broadcast over y)
        sPert.u = Am .* R.u;   % (Nmode x Ny x Nx) complex, physical peak field
        sPert.v = Am .* R.v;
        sPert.w = Am .* R.w;
        sPert.p = Am .* R.p;
        sPert.A = R.A;         % raw StabRes.A (= max|u'| peak); diagnostic
    else
        % No stored amplitude: assume u/v/w/p are already the physical fields.
        sPert.u = R.u;  sPert.v = R.v;  sPert.w = R.w;  sPert.p = R.p;
    end
    if isfield(R, 'omegavec'), sPert.omega = R.omegavec; end
    if isfield(R, 'betavec'),  sPert.beta  = R.betavec;  end
    if isfield(R, 'alpha'),    sPert.alpha = R.alpha;    end
    % reference scales (for normalizing plot axes): u_inf, length scale, Re
    if isfield(G, 'Uref'),     sPert.uref  = G.Uref;     end
    if isfield(G, 'lref'),     sPert.lref  = G.lref;     end
    if isfield(G, 'Re'),       sPert.Re    = G.Re;       end

    fprintf('importData: loaded %s — base flow %dx%d, %d perturbation mode(s)\n', ...
            inp.fieldsFile, Ny, Nx, size(sPert.u, 1));
end

% --- local helper: read midPlane.bin ---
% Binary layout written by the writeMidPlane function object:
%   header :  int32 nCols + int64 nRows           (12 bytes, native endian)
%   body   :  nRows rows × nCols doubles each     (row-major)
% Column meaning by nCols:
%   7  -> x, y, z, u, v, w, p             (DFP)
%   8  -> x, y, z, u, v, w, p, omz        (TTCP)
function T = readMidPlaneBinary(path)
    fid = fopen(path, 'r');
    if fid < 0
        error('importData:openFailed', 'Could not open %s', path);
    end
    cleaner = onCleanup(@() fclose(fid));
    nCols = fread(fid, 1, 'int32');
    nRows = fread(fid, 1, 'int64');
    data  = fread(fid, [double(nCols), double(nRows)], 'double')';
    if size(data,1) ~= nRows || size(data,2) ~= nCols
        error('importData:truncatedBin', ...
              'Truncated/short read in %s: header says %dx%d, got %dx%d.', ...
              path, nRows, nCols, size(data,1), size(data,2));
    end
    baseNames = {'x','y','z','u','v','w','p','omz'};
    if nCols < 7 || nCols > numel(baseNames)
        error('importData:badCols', ...
              'Unexpected nCols=%d in %s (expected 7 or 8).', nCols, path);
    end
    T = array2table(data, 'VariableNames', baseNames(1:nCols));
end

% --- local helper ---
function caseDir = caseDirFor(inp)
    switch inp.caseType
        case 'DFP'
            caseDir = fullfile(inp.airPowerRoot, 'PreProcessing', 'Modules', ...
                               'DirectFlatPlateModule');
        case 'TTCP'
            caseDir = fullfile(inp.airPowerRoot, 'PreProcessing', 'Modules', ...
                               'TunnelToCurvedPlateModule', 'AirfoilLECase');
        otherwise
            error('importData:badCaseType', ...
                  'Unknown caseType: ''%s''. Must be ''DFP'' or ''TTCP''.', ...
                  inp.caseType);
    end
end
