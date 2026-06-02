% inputs.m

% airPower root (auto-detected: assumes this file is at <airPower>/PostProcessing/inputs.m)
inp.airPowerRoot = fileparts(fileparts(mfilename('fullpath')));

% Which module's "post" outputs to import:
%   'DFP'  -> PreProcessing/Modules/DirectFlatPlateModule
%   'TTCP' -> PreProcessing/Modules/TunnelToCurvedPlateModule/AirfoilLECase
inp.caseType = 'DFP';
