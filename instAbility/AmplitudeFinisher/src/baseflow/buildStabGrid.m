function StabGrid = buildStabGrid(in)
% BUILDSTABGRID  OpenFOAM flat-plate base flow -> DeHNSSo StabGrid via gridgen.
%
%   StabGrid = buildStabGrid(in)
%
% Wraps DeHNSSo's main_gridgen: reads the OpenFOAM mid-plane csv
% (x,y,z,u,v,w,p,omz), resamples it onto the Malik-Chebyshev eta x equidistant
% xi grid (unstructured path), builds all metric terms and base-flow
% derivatives, and returns the StabGrid struct main() consumes. The result is
% cached to in.gridgen.outMat.
%
% in.gridgen.mode:
%   'auto'   run gridgen only if outMat is missing or older than the csv
%   'always' always re-run gridgen
%   'never'  require an existing outMat (error if absent), skip gridgen

gg = in.gridgen;

needBuild = strcmpi(gg.mode,'always');
if strcmpi(gg.mode,'auto')
    if ~isfile(gg.outMat)
        needBuild = true;
    else
        dcsv = dir(gg.baseFlowCsv);  dmat = dir(gg.outMat);
        if isempty(dcsv)
            error('buildStabGrid:noCsv','Base-flow csv not found: %s', gg.baseFlowCsv);
        end
        needBuild = dmat.datenum < dcsv.datenum;   % csv newer than cache
    end
elseif strcmpi(gg.mode,'never')
    if ~isfile(gg.outMat)
        error('buildStabGrid:noCache', ...
              'gridgen mode ''never'' but no StabGrid at %s', gg.outMat);
    end
    needBuild = false;
end

if needBuild
    if ~isfile(gg.baseFlowCsv)
        error('buildStabGrid:noCsv','Base-flow csv not found: %s', gg.baseFlowCsv);
    end
    fprintf('[buildStabGrid] Running DeHNSSo gridgen on %s ...\n', gg.baseFlowCsv);
    ggRoot = fullfile(in.dehnssoRoot, 'gridgen');
    addpath(ggRoot, genpath(fullfile(ggRoot, 'src')));

    [inFolder, inName, inExt] = fileparts(gg.baseFlowCsv);
    input = struct('folder', inFolder, 'filename', [inName inExt], ...
                   'format', 'csv', 'structured', false);
    [outFolder, outName, outExt] = fileparts(gg.outMat);
    output = struct('folder', outFolder, 'filename', [outName outExt]);

    StabGrid = main_gridgen(input, gg.params, output);
else
    fprintf('[buildStabGrid] Loading cached StabGrid: %s\n', gg.outMat);
    S = load(gg.outMat);
    if ~isfield(S,'StabGrid')
        error('buildStabGrid:badMat','%s has no StabGrid variable', gg.outMat);
    end
    StabGrid = S.StabGrid;
end
end
