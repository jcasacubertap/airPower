% inputs.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     PostProcessing — Inputs                          %
%                                                                      %
%   GENERAL INPUTS (top): which module to run + the shared data source.%
%   MODULE INPUTS (below): one block per module dispatched in run.m.   %
%   Only the block for the selected inp.task is used.                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% =====================================================================
%  GENERAL INPUTS
% =====================================================================

% airPower root (auto-detected: this file is at <airPower>/PostProcessing/inputs.m)
inp.airPowerRoot = fileparts(fileparts(mfilename('fullpath')));

% (1) Which post-processing module to run (dispatched in run.m):
%       'readData'           -> import the data only (base flow / perturbation),
%                               optionally overlay the validation plot
%       'productionAnalysis' -> Reynolds-Orr production of perturbation energy
inp.task = 'productionAnalysis';

% (2) Data source (shared by every module):
%   loadMode:
%     'loadBF'     -> base flow only, from a PreProcessing case's midPlane (inp.caseType)
%     'loadFields' -> base flow + perturbation, from io/fields/<inp.fieldsFile>
%                     (a DeHNSSo StabRes+StabGrid .mat)
inp.loadMode   = 'loadFields';
inp.caseType   = 'TTCP';              % (loadBF only)  'DFP' | 'TTCP'
inp.fieldsFile = 'dfp_linear.mat';    % (loadFields only) .mat inside io/fields

% =====================================================================
%  MODULE INPUTS   (only the selected task's block is used)
% =====================================================================

%% -- readData --
% Validation: when true, overlay the perturbation-amplitude comparison
% (max over y of |u'| vs x) against src/validation/<inp.validationFile>.
inp.validation     = false;
inp.validationFile = 'dataExtCrossA.mat';

%% -- productionAnalysis  (Reynolds-Orr) --
% Needs perturbation data (loadMode='loadFields'). No factor-2 correction: a
% stability code already combines the mode + its cc in the field array (see
% reynoldsOrrProdTermsFun header).
inp.ro.modeIdx   = [];        % [] -> all spanwise modes with beta>0 (MFD excluded)
inp.ro.xLim      = [];        % [xmin xmax] integration window (base-flow x units); [] -> full
inp.ro.yLim      = [];        % [ymin ymax] integration window; [] -> full
inp.ro.bufferFrac= 0.85;      % plot x-axis ends at this fraction of the domain (buffer start); 1 -> full
