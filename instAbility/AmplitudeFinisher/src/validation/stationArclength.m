function S_st = stationArclength(xc, in)
% STATIONARCLENGTH  Map PIV xc stations to arc-length from the leading edge.
%
%   S_st = stationArclength(xc, in)
%
% Reuses the DirectFlatPlate validation convention (dfpvalidation.jl,
% xc_to_xdfp): the PIV station label xc is a PERCENT-CHORD position, converted
% to a chordwise coordinate x_chord = (xc/100)*c and then to arc-length S along
% the surface by interpolating the reference boundary-layer table S(x) shipped
% in airfoilFlowData/<file>.mat (fields BL.S, BL.x, BL.c).
%
% Returns S_st in metres, in the SAME frame as the assembled base flow's BL.x
% (distance from the LE), so PIV and NPSE stations are directly comparable.
%
% The PIV station labels xc are interpreted as PERCENT-CHORD (x/c in %).

xc = double(xc(:)).';   % PIV station labels, percent-chord

% --- reference S(x) table -----------------------------------------------
refPath = fullfile(in.airPowerRoot, 'PreProcessing', 'io', 'airfoilFlowData', ...
                   in.validation.refFlowData);
if ~isfile(refPath)
    error('stationArclength:noRef', ...
          ['Reference flow-data .mat not found: %s\n(set ' ...
           'in.validation.refFlowData).'], refPath);
end
R  = load(refPath);
BL = R.BL;
xr = double(BL.x(:));           % chordwise coordinate [m]
Sr = double(BL.S(:));           % arc-length          [m]
c  = double(BL.c);              % chord               [m]

x_chord = (xc/100) * c;         % percent-chord -> chordwise coordinate [m]
S_st = interp1(xr, Sr, x_chord, 'linear', 'extrap');
end
