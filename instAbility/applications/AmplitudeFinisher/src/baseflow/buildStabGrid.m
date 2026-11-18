function [StabGrid, in] = buildStabGrid(in)
% BUILDSTABGRID  Load the DeHNSSo StabGrid (base flow) and set the frame offset.
%
%   [StabGrid, in] = buildStabGrid(in)
%
% The StabGrid is produced OUTSIDE this tool: run the DeHNSSo gridgen caller for
% your problem (in.caller), then paste the resulting DeHNSSo_input.mat into
% AmplitudeFinisher/io/ (in.stabGridFile). This function loads it and sets
% the streamwise frame offset in.match.xOffset so that
%     hns.x = StabGrid.x*lref + xOffset
% is arc-length S from the LE (same frame as the PIV stations):
%   * flat   -> xOffset = in.geom.xInlet (the flat-plate inlet distance from LE).
%   * curved -> xOffset = arc-length(LE -> grid origin), computed from the airfoil
%               geometry: invert x_phys(x/c) (as in PreProcessing airfoil_surface)
%               against the grid's Cartesian wall X, then map x/c -> S via the
%               reference BL table. This offset is constant to machine precision
%               because the grid wall-arc metric matches BL.S.
% There is deliberately NO CSV/gridgen/cache logic here.

    f = in.stabGridFile;
    if ~isfile(f)
        error('buildStabGrid:noFile', ...
              ['StabGrid not found:\n  %s\n\n' ...
               'Run the ''%s'' gridgen caller for the ''%s'' problem and paste its ' ...
               'DeHNSSo_input.mat there.'], f, in.caller, in.problem);
    end
    S = load(f);
    if ~isfield(S, 'StabGrid')
        error('buildStabGrid:noStabGrid', ...
              '%s has no ''StabGrid'' variable (is it a DeHNSSo gridgen output?).', f);
    end
    StabGrid = S.StabGrid;
    fprintf('[buildStabGrid] problem=%s  <-  %s\n', in.problem, f);

    % --- streamwise frame offset (per problem) ---------------------------
    switch lower(in.problem)
        case 'flat'
            in.match.xOffset = in.geom.xInlet;
        case 'curved'
            in.match.xOffset = local_curved_xoffset(StabGrid, in);
        otherwise
            error('buildStabGrid:problem', 'Unknown in.problem ''%s''.', in.problem);
    end
    fprintf('[buildStabGrid] xOffset = %.6g m  (%s frame)\n', in.match.xOffset, in.problem);
end

% -------------------------------------------------------------------------
function xOff = local_curved_xoffset(G, in)
% Arc-length from the LE to the curved grid's arc origin. Map each wall column
% to x/c by inverting the airfoil surface x_phys(x/c) (PreProcessing
% airfoil_surface), then to S via the reference BL(x -> S) table; the offset
% S - StabGrid.x*lref is constant, so return its mean.
    root = in.airPowerRoot;
    datF = fullfile(root, 'PreProcessing', 'io', 'airfoilGeometryData', in.geom.airfoilFile);
    if ~isfile(datF)
        error('buildStabGrid:noAirfoil', 'Airfoil geometry not found: %s', datF);
    end
    D = readmatrix(datF);
    xi = D(:,1); et = D(:,2);
    m = et >= 0; xi = xi(m); et = et(m);          % upper surface, LE -> TE
    [xi, p] = sort(xi); et = et(p);
    [xi, iu] = unique(xi); et = et(iu);
    c  = in.geom.chord;  al = in.geom.alphaDeg * pi/180;
    xf = linspace(0, 1, 4000).';  ef = interp1(xi, et, xf, 'spline');
    xphys = cos(al)*(xf - 0.5)*c + sin(al)*ef*c;  % airfoil_surface chordwise coord [m]
    mono  = xf <= 0.55;                            % monotonic front half (covers the window)

    if ~isfield(G, 'X')
        error('buildStabGrid:noCartesian', ...
              'Curved StabGrid lacks Cartesian ''X'' (needed for the x/c frame offset).');
    end
    [~, iw] = min(abs(G.y));                       % wall row (y = 0)
    Xw   = G.X(iw, :) * G.lref;                     % wall Cartesian x [m]
    xcCol = interp1(xphys(mono), xf(mono), Xw, 'linear', 'extrap');   % x/c fraction

    BLp = fullfile(root, 'PreProcessing', 'io', 'airfoilFlowData', in.validation.refFlowData);
    BL  = load(BLp);  BL = BL.BL;
    Sref = interp1(BL.x(:), BL.S(:), xcCol * c, 'linear', 'extrap');  % arc from LE [m]
    xOff = mean(Sref - G.x(:).' * G.lref);
end
