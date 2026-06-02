function [sBF, inp] = importBaseFlow(inp)
% importBaseFlow  Load a PreProcessing case's midPlane data and arrange it on a meshgrid.
%
%   [sBF, inp] = importBaseFlow(inp)
%
% Reads <case>/postProcessing/midPlane.{bin,csv} from the case selected
% by inp.caseType, and arranges its scalar fields onto a 2D meshgrid of
% (sBF.x, sBF.y). If both .bin and .csv are present, .bin wins.
%
% Returns sBF with fields (Ny x Nx matrices):
%   sBF.x, sBF.y           meshgrid coordinates
%   sBF.u, sBF.v, sBF.w    velocity components
%   sBF.p                  pressure
%   sBF.omz                z-vorticity (if present in source)
%
% Also sets inp.wallExtrap = true/false, depending on whether the imported
% midPlane file contained wall-extrapolation rows (u=v=w=0 added at wall
% face centres by the OpenFOAM writeMidPlane function object).
%
% Required inp fields:
%   inp.caseType       'DFP' or 'TTCP'
%   inp.airPowerRoot   absolute path to airPower

    caseDir = caseDirFor(inp);

    % --- read source: prefer binary, fall back to CSV ---
    csvPath = fullfile(caseDir, 'postProcessing', 'midPlane.csv');
    binPath = fullfile(caseDir, 'postProcessing', 'midPlane.bin');
    if isfile(binPath)
        T = readMidPlaneBinary(binPath);
        fprintf('importBaseFlow: read %s (%d rows)\n', binPath, height(T));
    elseif isfile(csvPath)
        T = readtable(csvPath);
        fprintf('importBaseFlow: read %s (%d rows)\n', csvPath, height(T));
    else
        error('importBaseFlow:notFound', ...
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
    fprintf('importBaseFlow: wallExtrap = %s (%d wall rows of %d total)\n', ...
            mat2str(inp.wallExtrap), nnz(isWall), height(T));

    % --- safety: the writeMidPlane function object must emit exactly one
    %     z-plane. If we see more than one distinct z value in the data,
    %     something on the OpenFOAM side has regressed (e.g. the z-tolerance
    %     widened, or the mesh got multiple cells per plane). Fail, so the
    %     user fixes the producer (system/controlDict, writeMidPlane).
    if ismember('z', T.Properties.VariableNames)
        zu = uniquetol(T.z, tol, 'DataScale', 1);
        if numel(zu) > 1
            error('importBaseFlow:multipleZPlanes', ...
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
        warning('importBaseFlow:nonRectangular', ...
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
        error('importBaseFlow:openFailed', 'Could not open %s', path);
    end
    cleaner = onCleanup(@() fclose(fid));
    nCols = fread(fid, 1, 'int32');
    nRows = fread(fid, 1, 'int64');
    data  = fread(fid, [double(nCols), double(nRows)], 'double')';
    if size(data,1) ~= nRows || size(data,2) ~= nCols
        error('importBaseFlow:truncatedBin', ...
              'Truncated/short read in %s: header says %dx%d, got %dx%d.', ...
              path, nRows, nCols, size(data,1), size(data,2));
    end
    baseNames = {'x','y','z','u','v','w','p','omz'};
    if nCols < 7 || nCols > numel(baseNames)
        error('importBaseFlow:badCols', ...
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
            error('importBaseFlow:badCaseType', ...
                  'Unknown caseType: ''%s''. Must be ''DFP'' or ''TTCP''.', ...
                  inp.caseType);
    end
end
