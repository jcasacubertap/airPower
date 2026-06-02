% run.m
% Entry point: load config, import post data into `sBF`.

inputs;
addpath(fullfile(fileparts(mfilename('fullpath')), 'Source'));
[sBF, inp] = importBaseFlow(inp);
