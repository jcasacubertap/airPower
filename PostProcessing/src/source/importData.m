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
%   sBF.u, sBF.v, sBF.w    velocity components
%   sBF.p                  pressure          (loadBF; absent in loadFields)
%   sBF.omz                z-vorticity       (loadBF, if present in source)
%   sBF.ux,uy,vx,vy,wx,wy  base-flow gradients — loadFields ONLY (read from the
%                          stability grid). loadBF provides none, for either mesh:
%                          both are multi-block structured and body-fitted, where a
%                          row step changes x and y together, so a Cartesian
%                          differentiateField call would be wrong. See
%                          assembleBodyFitted for what correct metrics would need.
%
% Coordinates come in one of two shapes, because the two sources hold the base
% flow on different kinds of grid:
%   sBF.x, sBF.y   rectilinear stability grid (loadFields), i.e. Cartesian x and y
%   sBF.s, sBF.n   wall-fitted pair on the body-fitted PreProcessing meshes
%                  (loadBF, both DFP and TTCP):
%                  s = arc length along the wall, 0 at the first exported station
%                  (the S* of the figures); n = wall-normal distance. Named (s,n)
%                  rather than (s,y) so that y NEVER means anything but Cartesian
%                  y, whichever grid is loaded.
%
% sBF.geom (loadBF on the body-fitted mesh) holds the physical Cartesian cell
% centres X, Y [m], co-registered with everything else. They are kept out of the
% top level because nothing plots against them — the figures all use (s,n) — and
% four coordinate matrices side by side only invites plotting the flow against
% the wrong pair. They are kept AT ALL because (X,Y) -> (s,n) is one-way: it
% discards where the wall is and how it is oriented, so the geometry cannot be
% rebuilt from (s,n). It is what the wall shape (including any modulation) lives
% in, what maps a station to x/c, and what the metric terms for curvilinear
% gradients are built from.
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

    % Mode selection is de-duplicated ONCE, here, because inp flows from this
    % single entry point to every plotter. A repeated index (inp.modeIdx =
    % [2 2], say) otherwise draws the same mode twice — overplotted curves and a
    % legend with two identical entries — in whichever plotter indexes modeIdx
    % directly. 'stable' keeps the user's ordering, so the first listed mode
    % still selects the fundamental in plotProfilesValidation.
    if isfield(inp, 'modeIdx') && ~isempty(inp.modeIdx)
        [u, keep] = unique(inp.modeIdx(:).', 'stable');
        if numel(u) < numel(inp.modeIdx)
            warning('importData:dupModeIdx', ...
                    'inp.modeIdx had %d repeated entr(y/ies); using [%s].', ...
                    numel(inp.modeIdx) - numel(u), num2str(u));
        end
        inp.modeIdx = u;  clear keep;
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
    % Nested geometry describes the same cells, so it flips with them or the two
    % stop being co-registered — the wall row would no longer line up with n = 0.
    if isfield(sBF, 'geom')
        gf = fieldnames(sBF.geom);
        for i = 1:numel(gf)
            if isequal(size(sBF.geom.(gf{i})), [Ny, Nx])
                sBF.geom.(gf{i}) = flipud(sBF.geom.(gf{i}));
            end
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

% --- loadBF: read a PreProcessing case's midPlane and put it back on a grid ---
% Reads <case>/postProcessing/midPlane.{bin,csv} from the case selected by
% inp.caseType. If both .bin and .csv are present, .bin wins. Also sets
% inp.wallExtrap depending on whether the file contained wall-extrapolation
% rows (u=v=w=0 added at wall face centres by the writeMidPlane function object).
%
% Both meshes are multi-block structured exports, so both go through
% assembleBodyFitted, which recovers the grid from the mesh ordering. Neither
% could be rebuilt by uniquetol on x and y: the cells are body-fitted, so every
% one carries its own (x,y) even though the (i,j) topology is perfect, and
% unique x/y sees ~N distinct values in each direction. The two differ only in
% how the blocks tile — TTCP one band of 6 equal-width blocks, DFP 2 bands of 6
% unequal ones — which assembleBodyFitted works out for itself.
% No path interpolates: every value returned is the cell value OpenFOAM wrote.
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

    % Absolute tolerance for the z-plane check below. The ASCII midPlane.csv
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

    sBF = assembleBodyFitted(T, isWall);
end

% --- DFP & TTCP: rebuild the multi-block body-fitted grid from the mesh ordering ---
% Returns, all Ny x Nx with row 1 = wall (orientFreestreamToWall then flips the
% whole set so row 1 is free-stream and row end the wall):
%   .s            streamwise arc length along the wall, 0 at the first exported
%                 station — the S* the figures are labelled with              [m]
%   .n            wall-normal distance measured up each grid column from the wall [m]
%   .geom.X/.Y    physical Cartesian cell centres [m] (same role as StabGrid.X/.Y),
%                 kept aside because no figure plots against them
%   .u .v .w .p .omz   the exported fields
%
% Base-flow gradients are deliberately NOT computed here. On a body-fitted grid
% a row step changes both x and y, so differentiateField(u, X, 2) returns
% u_xi / x_xi rather than du/dx; correct metric terms are a separate change.
% Nothing consumes loadBF gradients today (loadBF returns sPert = [], and run.m
% refuses reynoldsOrrProdTerms without perturbation data).
function sBF = assembleBodyFitted(T, isWall)

    Tw = T(isWall,  :);      % wall-face rows (u=v=w=0), one per station
    Ti = T(~isWall, :);      % interior cells
    if height(Ti) < 3
        error('importData:noInteriorCells', ...
              'midPlane holds %d interior cell(s) — nothing to assemble.', height(Ti));
    end

    % --- 1. cut the interior rows into streamwise sweeps ---
    % writeMidPlane appends cells in mesh order, which within a block means one
    % sweep along the block per wall-normal level. Neighbours inside a sweep sit
    % one streamwise spacing apart; the step from the end of a sweep back to the
    % start of the next spans the whole sweep, ~2 orders of magnitude more. That
    % jump marks the boundary, so no knowledge of the block layout is needed.
    step = hypot(diff(Ti.x), diff(Ti.y));
    mstep = median(step);
    if ~(mstep > 0)
        error('importData:degenerateOrder', ...
              'Cannot find streamwise sweeps: median cell-to-cell step is %g.', mstep);
    end
    brk  = find(step > 10 * mstep);
    sBeg = [1; brk + 1];
    sEnd = [brk; height(Ti)];
    sLen = sEnd - sBeg + 1;
    % Drop fragments at the exportHeight ceiling. The cut sits at a fixed
    % wall-normal distance, so a few columns can carry one extra cell right at it;
    % those arrive as a sweep far narrower than its neighbours. Compare LOCALLY,
    % not against a global width — block widths legitimately differ (DFP runs
    % 288/48/96/48/560/192), and only a fragment is narrow relative to both
    % neighbours. Keeping one would split a band and break the tiling.
    nb = [sLen(2:end); sLen(end)];  pv = [sLen(1); sLen(1:end-1)];
    frag = sLen < 0.5 * min(nb, pv);
    if any(frag)
        fprintf(['importData: dropped %d cell(s) in %d partial sweep(s) at the ' ...
                 'exportHeight ceiling (grid kept rectangular)\n'], ...
                sum(sLen(frag)), nnz(frag));
        sBeg = sBeg(~frag);  sEnd = sEnd(~frag);  sLen = sLen(~frag);
    end

    nSweep = numel(sBeg);
    if nSweep < 3
        error('importData:tooFewSweeps', ...
              ['Found only %d sweep(s) in the midPlane export — the mesh ordering ' ...
               'is not what assembleBodyFitted expects. Check writeMidPlane in ' ...
               '<case>/system/controlDict.'], nSweep);
    end

    % --- 2. group the sweeps into mesh blocks ---
    % blockMesh numbers cells block by block, so one block's sweeps are contiguous,
    % all the same width, and their start points climb slowly away from the wall.
    % A boundary is therefore a change of width OR a jump in where the sweep starts
    % (the next block sits somewhere else entirely). Both tests are needed: TTCP's
    % blocks are all the same width and are caught by the jump, while DFP's differ
    % in width (288/48/96/48/560/192) and some are caught by the width change.
    p0x = Ti.x(sBeg);  p0y = Ti.y(sBeg);
    d0  = hypot(diff(p0x), diff(p0y));
    isBnd = (diff(sLen) ~= 0) | (d0 > 10 * median(d0));
    bBeg  = [1; find(isBnd) + 1];
    bEnd  = [bBeg(2:end) - 1; nSweep];
    nBlk  = numel(bBeg);

    W = sLen(bBeg);                     % stations per block
    L = bEnd - bBeg + 1;                % wall-normal levels per block

    % --- 3. place each block in the global grid ---
    % Blocks tile a rectangle, but in BOTH directions: DFP is 6 blocks across x by
    % 2 bands in y ($NyB / $NyT in its blockMeshDict), TTCP a single band of 6.
    % A band is the set of blocks spanning the same wall-normal interval, so group
    % by that, order the bands wall-outwards, and order each band streamwise.
    % A band is a run of CONSECUTIVE blocks carrying the same number of levels:
    % blockMesh numbers a whole band before starting the next ($NyB blocks 0-5,
    % then $NyT blocks 6-11 in the DFP dict), so the run structure identifies them
    % directly. Matching wall-normal extents numerically does not work — the bump
    % deforms the mesh, so blocks in one band no longer share an extent.
    band  = cumsum([1; L(2:end) ~= L(1:end-1)]);
    nBand = band(end);

    dLo = zeros(nBlk,1);
    for b = 1:nBlk
        lo = sBeg(bBeg(b)) : sEnd(bBeg(b));          % sweep nearest the wall
        dLo(b) = mean(wallDist(Ti.x(lo), Ti.y(lo), Tw));
    end
    [~, bandOrder] = sort(arrayfun(@(k) min(dLo(band == k)), 1:nBand));   % wall outwards

    Nx = sum(W(band == bandOrder(1)));
    Ny = 0;
    for k = bandOrder
        if sum(W(band == k)) ~= Nx
            error('importData:blocksNotRectangular', ...
                  ['Mesh blocks do not tile a rectangle: band %d spans %d stations ' ...
                   'but the wall band spans %d.'], k, sum(W(band == k)), Nx);
        end
        Ny = Ny + L(find(band == k, 1));
    end
    fprintf('importData: %d blocks in %d band(s) -> %dx%d grid\n', nBlk, nBand, Ny, Nx);

    % --- 4. reshape each block and lay it into place ---
    % Inside a block the rows run sweep by sweep, so an (W x L) reshape transposed
    % gives level-by-station directly.
    names  = {'x','y','u','v','w','p','omz'};
    names  = names(ismember(names, Ti.Properties.VariableNames));
    G = struct();
    for k = 1:numel(names); G.(names{k}) = nan(Ny, Nx); end

    r0 = 0;
    for k = bandOrder
        % Blocks keep the mesh's own order within a band. Sorting them by mean x
        % looks tempting but is wrong: it reorders blocks without reversing the
        % columns INSIDE each one, and TTCP's blocks run x-descending internally,
        % so the result is a sawtooth. Step 5 checks the assembled wall instead.
        ib = find(band == k);
        c0 = 0;
        for b = ib(:).'
            rows = sBeg(bBeg(b)) : sEnd(bEnd(b));
            for m = 1:numel(names)
                G.(names{m})(r0 + (1:L(b)), c0 + (1:W(b))) = ...
                    reshape(Ti.(names{m})(rows), W(b), L(b)).';
            end
            c0 = c0 + W(b);
        end
        r0 = r0 + L(ib(1));
    end

    % --- 5. attach the wall row, and use it to check the station order ---
    % wallExtrapolation writes one row per wall face in wall-face order, which is
    % the order the blocks tile. Each wall point must therefore sit directly
    % below the matching column's first cell; if it does not, the blocks came out
    % in an order this function does not understand and the grid would be
    % scrambled, so fail rather than return a plausible-looking wrong field.
    haveWall = height(Tw) == Nx;
    if height(Tw) > 0 && ~haveWall
        warning('importData:wallRowCount', ...
                ['midPlane holds %d wall row(s) but the grid has %d station(s); ' ...
                 'the wall row is dropped and wall-normal distance is measured ' ...
                 'from the first cell instead.'], height(Tw), Nx);
    end
    if haveWall
        % Match wall faces to columns rather than assuming a shared order: the
        % blocks are laid out streamwise here, which need not be the order
        % writeMidPlane emitted the wall faces in. A correct grid gives a
        % bijection with every wall point within a cell or so of its column's
        % first cell; a scrambled one gives neither, so both are checked.
        [off, pick] = min(hypot(Tw.x(:) - G.x(1,:), Tw.y(:) - G.y(1,:)), [], 1);
        if numel(unique(pick)) ~= Nx || max(off) > 5 * mstep
            error('importData:stationOrderMismatch', ...
                  ['Wall rows do not map one-to-one onto the reconstructed ' ...
                   'stations (%d distinct of %d; worst offset %.3g m against a ' ...
                   'streamwise spacing of %.3g m). The block layout is not what ' ...
                   'assembleBodyFitted expects.'], ...
                  numel(unique(pick)), Nx, max(off), mstep);
        end
        for k = 1:numel(names)
            f = names{k};
            G.(f) = [Tw.(f)(pick).'; G.(f)];   % wall becomes row 1, in column order
        end
        Ny = Ny + 1;

        % The assembled wall must be a continuous path: neighbouring columns one
        % streamwise spacing apart. Blocks laid out in the wrong order still pass
        % the bijection above (every column keeps a unique nearest wall face) but
        % leave a jump at each block seam, which would silently inflate the arc
        % length s computed from this row.
        wstep = hypot(diff(G.x(1,:)), diff(G.y(1,:)));
        if max(wstep) > 5 * median(wstep)
            error('importData:wallNotContiguous', ...
                  ['The assembled wall jumps %.3g m between stations %d and %d ' ...
                   '(typical spacing %.3g m), so the blocks are not laid out ' ...
                   'head-to-tail along the wall.'], ...
                  max(wstep), find(wstep == max(wstep), 1), ...
                  find(wstep == max(wstep), 1) + 1, median(wstep));
        end
    end

    % --- 5b. orient the stations inflow -> outflow ---
    % The block order is the mesh's, and it need not follow the flow: on the TTCP
    % C-grid the wall faces run from xi/c ~ 0.5 back to ~ 0.02, i.e. AGAINST it.
    % Decide physically, exactly as the wall-normal order is decided downstream —
    % project the free-stream velocity (row end here) onto the wall tangent that
    % points from column i to i+1 (row 1 is the wall). If the flow opposes it,
    % reverse every column so column 1 is the most upstream station. s then starts
    % at the numerical inflow as advertised, the boundary layer thickens with s
    % instead of thinning, and a downstream-fraction cut (inp.plot.bufferFrac)
    % trims the outflow buffer rather than the inlet.
    tx = gradient(G.x(1,:));  ty = gradient(G.y(1,:));
    tl = hypot(tx, ty);  tx = tx ./ tl;  ty = ty ./ tl;
    proj = G.u(end,:) .* tx + G.v(end,:) .* ty;
    if mean(proj, 'omitnan') < 0
        for k = 1:numel(names)
            G.(names{k}) = fliplr(G.(names{k}));
        end
        fprintf(['importData: reversed station order (mesh ran against the flow) ' ...
                 '-> column 1 is the inflow\n']);
    else
        fprintf('importData: stations already run inflow -> outflow\n');
    end

    % --- 6. coordinates ---
    % Wall-normal distance up each column, 0 at the wall (row 1 here), and
    % streamwise arc length along that wall row, 0 at the first exported station.
    % Same construction as plotCoords, so the two agree exactly. Without a wall
    % row, n = 0 falls on the first interior cell instead. Note s is the WALL
    % station broadcast up the column, not the local arc length at height n.
    seg   = hypot(diff(G.x, 1, 1), diff(G.y, 1, 1));
    sBF.n = [zeros(1, Nx); cumsum(seg, 1)];
    sw    = [0, cumsum(hypot(diff(G.x(1,:)), diff(G.y(1,:))))];
    sBF.s = repmat(sw, Ny, 1);

    % Physical Cartesian cell centres, one level down: the figures never use them,
    % but they are the only record of where the wall actually is (see the header).
    sBF.geom.X = G.x;
    sBF.geom.Y = G.y;

    for k = 1:numel(names)
        f = names{k};
        if ~ismember(f, {'x','y'}); sBF.(f) = G.(f); end
    end
end

% --- local helper: distance from each point to the nearest wall face ---
% Used only to order and group the mesh blocks, so a nearest-centre distance is
% accurate enough even where the wall curves. Without wall rows, fall back to y.
function d = wallDist(x, y, Tw)
    if isempty(Tw)
        d = y(:);  return;
    end
    d = min(hypot(x(:) - Tw.x(:).', y(:) - Tw.y(:).'), [], 2);
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

% --- local helper: PreProcessing case directory for the selected caseType ---
% Segments are spelled exactly as they are on disk, and resolveCase enforces
% that. macOS resolves paths case-INSENSITIVELY and Linux does not, so a folder
% rename that is not mirrored here keeps every macOS run working while breaking
% every Linux run — which is exactly what happened when the tree moved to
% lower-case names and these four segments were left behind.
function caseDir = caseDirFor(inp)
    switch inp.caseType
        case 'DFP'
            segs = {'PreProcessing', 'modules', 'directFlatPlateModule'};
        case 'TTCP'
            segs = {'PreProcessing', 'modules', 'tunnelToCurvedPlateModule', ...
                    'airfoilLECase'};
        otherwise
            error('importData:badCaseType', ...
                  'Unknown caseType: ''%s''. Must be ''DFP'' or ''TTCP''.', ...
                  inp.caseType);
    end
    caseDir = resolveCase(inp.airPowerRoot, segs);
end

% --- walk a path one segment at a time, checking each against the real listing ---
% An exact match is taken silently. A segment that differs only in case still
% resolves — so a rename never hard-stops the user — but warns, on BOTH
% platforms, naming the spelling to correct here. A missing segment is an error.
function p = resolveCase(root, segs)
    p = root;
    for k = 1:numel(segs)
        d = dir(p);
        names = {d([d.isdir]).name};
        if any(strcmp(names, segs{k}))
            p = fullfile(p, segs{k});
            continue;
        end
        hit = names(strcmpi(names, segs{k}));
        if isempty(hit)
            error('importData:caseDirNotFound', ...
                  'No directory named ''%s'' under %s.', segs{k}, p);
        end
        warning('importData:caseDirSpelling', ...
                ['Directory ''%s'' is spelled ''%s'' on disk. Resolved, but ' ...
                 'this only works on a case-insensitive filesystem — update ' ...
                 'caseDirFor in importData.m so Linux works too.'], ...
                segs{k}, hit{1});
        p = fullfile(p, hit{1});
    end
end
