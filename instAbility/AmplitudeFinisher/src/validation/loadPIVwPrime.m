function piv = loadPIVwPrime(valRoot, gen, caseId, in)
% LOADPIVWPRIME  Load the experimental PIV w' target for a Gen/Case.
%
%   piv = loadPIVwPrime(valRoot, gen, caseId, in)
%
% Reads Validation/Gen{gen}/Experimental/Case{caseId}/*.mat. The file holds a
% struct (in.validation.structVar, default 'output') of 1 x nStations cell
% arrays. For each spanwise-Fourier mode k = 1..nModes it extracts the
% wall-normal RMS profile of w': w_pert_m_prof_rms_0k(y_prof_rms_0k), at the
% streamwise stations xc.
%
% Coordinates/velocities are converted to SI (in.validation.*Scale). Because
% each mode is sampled on its own y-grid, the per-mode profiles are
% interpolated onto a common per-station y-grid (mode 1's) so they can be
% combined; the combined RMS profile (Parseval, sqrt(sum_k w_k^2)) is the
% optimization target, the per-mode profiles are kept for diagnostics.
%
% Returns a struct:
%   piv.x        1 x ns          streamwise stations                 [m]
%   piv.y        1 x ns cell     wall-normal coordinate per station  [m]
%   piv.rms      1 x ns cell     combined RMS w' profile (target)    [m/s]
%   piv.mode     1 x ns cell     [ny x nModes] per-mode w' profiles  [m/s]
%   piv.file     source .mat path
%   piv.nModes   number of modes loaded

caseDir = fullfile(valRoot, sprintf('Gen%d', gen), 'Experimental', ...
                   sprintf('Case%d', caseId));
if ~isfolder(caseDir)
    error('loadPIVwPrime:noCase', 'Case dir not found: %s', caseDir);
end

% --- locate & load the .mat ---------------------------------------------
if ~isempty(in.validation.matFile)
    matPath = fullfile(caseDir, in.validation.matFile);
else
    hits = dir(fullfile(caseDir, '*.mat'));
    if isempty(hits)
        error('loadPIVwPrime:noMat', 'No .mat files in %s', caseDir);
    end
    matPath = fullfile(caseDir, hits(1).name);
end
S = load(matPath);

sv = in.validation.structVar;
if ~isfield(S, sv)
    error('loadPIVwPrime:noStruct', ...
          '%s has no struct ''%s''. Found: %s', ...
          matPath, sv, strjoin(fieldnames(S), ', '));
end
o = S.(sv);

local_assert_fields(o, {in.validation.xField});

% How many spanwise-Fourier modes the file actually contains: count the
% leading per-mode fields (w_pert_m_prof_rms_0k / y_prof_rms_0k) present.
fn = fieldnames(o);
nModes = 0;
while nModes < numel(in.validation.wFields) ...
        && ismember(in.validation.wFields{nModes+1}, fn) ...
        && ismember(in.validation.yFields{nModes+1}, fn)
    nModes = nModes + 1;
end
if nModes == 0
    error('loadPIVwPrime:noModes', ...
          'No per-mode PIV w'' fields found (looked for %s ...). Available: %s', ...
          in.validation.wFields{1}, strjoin(fn, ', '));
end

% --- per-station assembly ------------------------------------------------
ns  = numel(o.(in.validation.xField));
Lsc = in.validation.lengthScale;
Vsc = in.validation.velScale;

piv.file   = matPath;
piv.nModes = nModes;
piv.x      = zeros(1, ns);
piv.y      = cell(1, ns);
piv.rms    = cell(1, ns);       % sqrt(sum_k mode_k^2)  (fallback target)
piv.mode   = cell(1, ns);       % per-mode RMS profiles (diagnostics)
piv.rmsFull = cell(1, ns);      % raw full-plane w' RMS  (peak-match target)

hasTot = isfield(in.validation,'totField') && ~isempty(in.validation.totField) ...
         && isfield(o, in.validation.totField);

% Streamwise stations -> arc-length from LE (same frame as base-flow BL.x),
% reusing the DFP percent-chord -> S(x) mapping (stationArclength).
xc_raw = zeros(1, ns);
for i = 1:ns; xc_raw(i) = double(o.(in.validation.xField){i}); end
piv.x  = stationArclength(xc_raw, in);   % arc-length S from LE [m]
piv.xc = xc_raw;                         % raw station labels (percent chord)

for i = 1:ns
    % Common y-grid for this station: mode 1's wall-normal coordinate.
    y1 = double(o.(in.validation.yFields{1}){i}(:)) * Lsc;
    ny = numel(y1);
    modeMat = zeros(ny, nModes);

    for k = 1:nModes
        yk = double(o.(in.validation.yFields{k}){i}(:)) * Lsc;
        wk = double(o.(in.validation.wFields{k}){i}(:)) * Vsc;
        if isequal(yk, y1)
            modeMat(:, k) = wk;
        else
            % Different per-mode grids -> interpolate onto y1.
            modeMat(:, k) = interp1(yk, wk, y1, 'linear', 'extrap');
        end
    end

    piv.y{i}    = y1;
    piv.mode{i} = modeMat;
    % Per-mode combined RMS (Parseval): sqrt(sum_k w_k^2).
    piv.rms{i}  = sqrt(sum(modeMat.^2, 2));

    % Raw full-plane w' RMS: spanwise (z) RMS of the mean-removed total plane.
    % z varies along dim 2; the plane rows align with the y1 grid.
    if hasTot
        wtot  = double(o.(in.validation.totField){i}) * Vsc;   % ny x nz
        wfluc = wtot - mean(wtot, 2);                          % remove spanwise mean
        rf    = sqrt(mean(wfluc.^2, 2));                       % ny x 1 RMS over z
        if numel(rf) == numel(y1)
            piv.rmsFull{i} = rf;
        else
            piv.rmsFull{i} = interp1(linspace(0,1,numel(rf)), rf, ...
                                     linspace(0,1,numel(y1)), 'linear', 'extrap').';
        end
    else
        piv.rmsFull{i} = piv.rms{i};   % fall back to mode-sum if no raw plane
    end
end
end

% -------------------------------------------------------------------------
function local_assert_fields(o, names)
missing = names(~ismember(names, fieldnames(o)));
if ~isempty(missing)
    error('loadPIVwPrime:missingVar', ...
          'PIV struct missing fields: %s. Available: %s', ...
          strjoin(missing, ', '), strjoin(fieldnames(o), ', '));
end
end
