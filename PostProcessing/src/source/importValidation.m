function sVal = importValidation(inp)
% importValidation  Load reference perturbation-amplitude curves for validation.
%
%   sVal = importValidation(inp)
%
% Reads src/validation/<inp.validationFile> (a MATLAB v7.3 / HDF5 .mat that
% holds a struct 'g') and returns the chordwise amplitude curves used to
% validate the loaded sPert. Only the small amplitude datasets are read (via
% h5read) — not the full 3D fields — so the (large) file is not fully loaded.
%
% Returns sVal with (modes along rows, to match sPert):
%   sVal.x                          (1 x Nxr) reference chordwise coord (g.Xun)
%   sVal.AMaxU/AMaxV/AMaxW/AMaxP    (Nmode x Nxr) max over y of |.'| per mode
%   sVal.lref, sVal.uref, sVal.Re   reference scales / Reynolds number (g.Reg)
%
% Required inp fields:
%   inp.airPowerRoot     absolute path to airPower
%   inp.validationFile   name of the .mat inside src/validation

    if ~isfield(inp, 'validationFile') || isempty(inp.validationFile)
        error('importValidation:noFile', ...
              'inp.validationFile is required to load validation data.');
    end

    matPath = fullfile(inp.airPowerRoot, 'PostProcessing', 'src', 'validation', ...
                       inp.validationFile);
    if ~isfile(matPath)
        error('importValidation:notFound', 'Validation file not found: %s', matPath);
    end

    % NOTE: h5read returns MATLAB-v7.3 arrays transposed (dims reversed), so the
    % AMax* curves come back as (Nxr x Nmode) and are transposed to (Nmode x Nxr).
    sVal.x     = h5read(matPath, '/g/Xun');   sVal.x = sVal.x(:).';   % 1 x Nxr
    sVal.AMaxU = h5read(matPath, '/g/AMaxU').';
    sVal.AMaxV = h5read(matPath, '/g/AMaxV').';
    sVal.AMaxW = h5read(matPath, '/g/AMaxW').';
    sVal.AMaxP = h5read(matPath, '/g/AMaxP').';
    sVal.lref  = h5read(matPath, '/g/lref');
    sVal.uref  = h5read(matPath, '/g/uref');
    sVal.Re    = h5read(matPath, '/g/Reg');   % reference Reynolds number (g.Reg)

    fprintf('importValidation: loaded %s — %d reference mode(s), %d chordwise points\n', ...
            inp.validationFile, size(sVal.AMaxU, 1), numel(sVal.x));
end
