% run.m
% Entry point: load config, import post data into `sBF` (and `sPert` for loadFields).

clear; clc;

% Load configuration settings
inputs;

% Add paths
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'src')));

% Import data
[sBF, sPert, inp] = importData(inp);

% Plot
if isfield(inp, 'validation') && inp.validation
    sVal = importValidation(inp);
    plotPerturbationAmplitude(sBF, sPert, inp, sVal);
end
