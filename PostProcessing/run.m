% run.m
% Entry point: load config, import the data, run the selected post-processing
% module (inp.task). Data import is shared by every module.

clear; clc;

% Load configuration
inputs;

% Add paths
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'src')));

% Import data (base flow, and perturbation for loadFields) — shared by all modules
[sBF, sPert, inp] = importData(inp);

% Dispatch on the selected module
switch lower(inp.task)

    case 'readdata'
        % Data only. Optionally overlay the perturbation-amplitude validation.
        if isfield(inp, 'validation') && inp.validation
            sVal = importValidation(inp);
            plotPerturbationAmplitude(sBF, sPert, inp, sVal);
        end

    case 'productionanalysis'
        % Reynolds-Orr production of perturbation kinetic energy.
        if isempty(sPert)
            error('run:noPert', ...
                  ['productionAnalysis needs perturbation data — ', ...
                   'set inp.loadMode = ''loadFields''.']);
        end
        sRO = reynoldsOrrProdTermsFun(sBF, sPert, inp);
        plotReynoldsOrrProd(sRO, sBF, inp);     % total P    -> io/plotting/reynoldsOrrProd.png
        plotReynoldsOrrDecomp(sRO, sBF, inp);   % I1..I4     -> io/plotting/reynoldsOrrDecomp.png

    otherwise
        error('run:task', ...
              'Unknown inp.task ''%s'' (expected ''readData'' | ''productionAnalysis'').', ...
              inp.task);
end
