% inputs.m

% airPower root (auto-detected: assumes this file is at <airPower>/PostProcessing/inputs.m)
inp.airPowerRoot = fileparts(fileparts(mfilename('fullpath')));

% How to load the data:
%   'loadBF'     -> base flow only, from PreProcessing midPlane (uses inp.caseType)
%   'loadFields' -> base flow + perturbation together, from io/fields/<inp.fieldsFile>.
%                   Here the base flow has been interpolated onto the stability
%                   grid, so it differs from the midPlane base flow.
inp.loadMode = 'loadFields';

% Which module's "post" outputs to import (loadBF only):
%   'DFP'  -> PreProcessing/Modules/DirectFlatPlateModule
%   'TTCP' -> PreProcessing/Modules/TunnelToCurvedPlateModule/AirfoilLECase
inp.caseType = 'DFP';

% Name of the .mat file inside io/fields to load (loadFields only):
inp.fieldsFile = 'Casacuberta2022.mat';

% Validation: when true, plot the perturbation-amplitude comparison
% (max over y of |u'| vs x) for the loaded sPert. The reference curves
% from src/validation/<inp.validationFile> will be overlaid later.
inp.validation     = true;
inp.validationFile = 'dataExtCrossA.mat';
