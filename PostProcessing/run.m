% run.m
% Entry point: load config, then either COMPUTE the selected post-processing
% analysis (inp.task) or — with inp.loadAnalysis=true — LOAD a previously-saved
% io/output bundle and re-plot without recomputing.

clear; clc;

% Load configuration — generated from airPower/inputs.jl. This entry point works
% BOTH ways, and always honours inputs.jl:
%   (1) via the dispatcher  — `julia run.jl PostProcessing <task>` writes a fresh
%       inputs_gen.m (with the CLI task) and sets AIRPOWER_PP_DISPATCH before
%       launching MATLAB here. We then load that config as-is (CLI task wins).
%   (2) run DIRECTLY in MATLAB — to load/inspect data interactively. No dispatch
%       flag is set, so we (re)generate inputs_gen.m from inputs.jl first, so a
%       direct run always mirrors the current inputs.jl (task = its default).
here     = fileparts(mfilename('fullpath'));
genFile  = fullfile(here, 'inputs_gen.m');

if ~isempty(getenv('AIRPOWER_PP_DISPATCH'))
    % (1) dispatcher run — config already written with the CLI task.
    if ~isfile(genFile)
        error('run:noConfig', 'inputs_gen.m not found in %s (dispatcher did not write it).', here);
    end
else
    % (2) direct MATLAB run — re-sync inputs_gen.m from inputs.jl via the dispatcher.
    fprintf('run: syncing inputs_gen.m from inputs.jl ...\n');
    % Locate julia: AIRPOWER_JULIA env → julia on PATH → ~/.juliaup/bin/julia.
    jl = getenv('AIRPOWER_JULIA');
    if isempty(jl) || ~isfile(jl)
        [s0, p0] = system('command -v julia');
        if s0 == 0 && ~isempty(strtrim(p0))
            jl = strtrim(p0);
        else
            cand = fullfile(getenv('HOME'), '.juliaup', 'bin', 'julia');
            if isfile(cand)
                jl = cand;
            else
                error('run:noJulia', ...
                      ['julia not found (needed to sync inputs_gen.m from inputs.jl).\n' ...
                       'Set AIRPOWER_JULIA to the julia binary, or run via the dispatcher.']);
            end
        end
    end
    cmd = sprintf('"%s" "%s" PostProcessing config', jl, fullfile(here, '..', 'run.jl'));
    [st, out] = system(cmd);
    if st ~= 0 || ~isfile(genFile)
        error('run:genConfig', ...
              ['Could not sync inputs_gen.m from inputs.jl.\n' ...
               'Fix julia (set AIRPOWER_JULIA), or run via the dispatcher:\n' ...
               '  julia run.jl PostProcessing <importData|reynoldsOrrProdTerms>\n%s'], out);
    end
end
addpath(here);
inputs_gen;   % sets `inp` (from the PostProcessing block of airPower/inputs.jl)

% Add paths
addpath(genpath(fullfile(here, 'src')));

if isfield(inp, 'loadAnalysis') && inp.loadAnalysis
    % ===== LOAD an existing analysis bundle from io/output and re-plot =====
    [~, stabname] = fileparts(inp.fieldsFile);
    bundle = fullfile(here, 'io', 'output', [stabname '_' inp.task '.mat']);
    if ~isfile(bundle)
        error('run:noBundle', ...
              ['loadAnalysis=true but bundle not found:\n  %s\n' ...
               'Run with loadAnalysis=false first to compute it.'], bundle);
    end
    S = load(bundle);
    fprintf('loaded analysis bundle <- %s\n', bundle);
    % Expose the bundle contents under the same names as a compute run, so the
    % workspace is identical either way. Keep the current `inp` (may carry updated
    % plot options); the bundle's own inp stays inside S.
    sBF = S.sBF; sPert = S.sPert; sRO = S.sRO;
    switch lower(inp.task)
        case {'productionanalysis', 'reynoldsorrprodterms'}
            plotReynoldsOrrProd(sRO, sBF, inp);
            plotReynoldsOrrDecomp(sRO, sBF, inp);
        otherwise
            error('run:loadTask', ...
                  'loadAnalysis only supports reynoldsOrrProdTerms (got ''%s'').', inp.task);
    end

else
    % ===== COMPUTE: import data, run the module, save + plot =====
    [sBF, sPert, inp] = importData(inp);

    switch lower(inp.task)

        case {'readdata', 'importdata'}
            % Base-flow plots always; perturbation plots too when loadFields
            % provided sPert. Figures are saved as PNGs only in the terminal /
            % dispatcher workflow (AIRPOWER_PP_DISPATCH set) — under
            % io/plotting/<caseFolder>; a direct MATLAB run leaves them open.
            if ~isempty(getenv('AIRPOWER_PP_DISPATCH'))
                switch upper(inp.caseType)
                    case 'DFP';  caseFolder = 'directFlatPlate';
                    case 'TTCP'; caseFolder = 'tunnelToCurvedPlate';
                    otherwise;   caseFolder = 'directFlatPlate';
                end
                savedir = fullfile(inp.airPowerRoot, 'PostProcessing', 'io', 'plotting', caseFolder);
                if ~exist(savedir, 'dir'); mkdir(savedir); end
            else
                savedir = '';   % direct MATLAB: show figures, do not save
            end

            plotBaseFlow(sBF, inp, savedir);                       % base flow (both modes)
            if ~isempty(sPert)
                plotPerturbationFields(sPert, sBF, inp, savedir);  % perturbation shapes + amplitude
                % wall-normal profiles at stations: u (always, non-dim); and,
                % when inp.validation, the dimensional w comparison vs PIV.
                plotProfiles(sBF, sPert, inp, savedir);
            end

        case {'productionanalysis', 'reynoldsorrprodterms'}
            % Reynolds-Orr production of perturbation kinetic energy.
            if isempty(sPert)
                error('run:noPert', ...
                      ['productionAnalysis needs perturbation data — ', ...
                       'set inp.loadMode = ''loadFields''.']);
            end
            sRO = reynoldsOrrProdTerms(sBF, sPert, inp);
            plotReynoldsOrrProd(sRO, sBF, inp);     % total P    -> io/plotting/reynoldsOrrProd.png
            plotReynoldsOrrDecomp(sRO, sBF, inp);   % I1..I4     -> io/plotting/reynoldsOrrDecomp.png

            % Self-contained analysis bundle: production + base flow + perturbation
            % (+ inp as metadata), keyed to the Stab run it came from.
            outdir = fullfile(inp.airPowerRoot, 'PostProcessing', 'io', 'output');
            if ~exist(outdir, 'dir'); mkdir(outdir); end
            [~, stabname] = fileparts(inp.fieldsFile);
            outfile = fullfile(outdir, [stabname '_' inp.task '.mat']);
            save(outfile, 'sRO', 'sBF', 'sPert', 'inp');
            fprintf('saved analysis bundle -> %s\n', outfile);

        otherwise
            error('run:task', ...
                  'Unknown inp.task ''%s'' (expected importData | reynoldsOrrProdTerms).', ...
                  inp.task);
    end
end

% Tidy the workspace: drop the scratch (paths, handles, loop temporaries) and keep
% only the analysis outputs, so a direct MATLAB run leaves a clean workspace.
clearvars -except sBF sPert sRO inp;
