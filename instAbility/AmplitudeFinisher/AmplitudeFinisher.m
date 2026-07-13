%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       AmplitudeFinisher                               %
%                                                                       %
%   Find the disturbance amplitude (Stab.A0_fund) that makes the        %
%   DeHNSSo Harmonic Navier-Stokes solution's w'-perturbation match     %
%   the PIV, for an OpenFOAM flat-plate base flow.                      %
%                                                                       %
%   For a new base flow: set the MAIN INPUTS in                         %
%   inputsAmplitudeFinisher.m, then just run this file.                 %
%                                                                       %
%   Pipeline:                                                           %
%     base-flow CSV --gridgen--> StabGrid                               %
%     PIV Gen/Case  ----------->  w' target (window x/c)                %
%     -> findAmplitude (linear seed -> nonlinear peak-match sweep)      %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; format long

%% Paths
here = fileparts(mfilename('fullpath'));
addpath(here, genpath(fullfile(here, 'src')));
AIRPOWER = fullfile(here, '..', '..');               % airPower root
in = inputsAmplitudeFinisher(AIRPOWER);
addpath(genpath(fullfile(in.dehnssoRoot, 'src')));   % DeHNSSo solver

%% 1. Base flow -> StabGrid
fprintf('\n[1/3] Building StabGrid from base flow ...\n');
StabGrid = buildStabGrid(in);

%% 2. PIV target for the requested Gen/Case
if isempty(in.VAL.Gen) || isempty(in.VAL.Case)
    [gen, caseId] = resolveGenCase(in.airPowerInputs);
else
    gen = in.VAL.Gen; caseId = in.VAL.Case;
end
fprintf('[2/3] Loading PIV: Gen%d / Case%d ...\n', gen, caseId);
piv = loadPIVwPrime(in.validation.rootDir, gen, caseId, in);
in.validation.nModes = piv.nModes;   % HNS diagnostics match the PIV modes present

% Select the stations in the chordwise window and truncate the domain to them.
sel = find(piv.xc >= in.match.xcWindow(1) & piv.xc <= in.match.xcWindow(2));
if isempty(sel)
    error('AmplitudeFinisher:emptyWindow', ...
          'No PIV stations with x/c in [%g %g]. Available x/c: %s', ...
          in.match.xcWindow(1), in.match.xcWindow(2), num2str(piv.xc));
end
in.match.stationIdx = sel;
xEnd = piv.x(sel(end)) / in.match.bufferFrac;
[StabGrid, tinfo] = truncateStabGrid(StabGrid, xEnd, in);
fprintf('      Matching x/c = %s ; domain nx %d -> %d\n', ...
        num2str(piv.xc(sel)), tinfo.nxOld, tinfo.nxNew);

%% 3. Find the matching amplitude
fprintf('[3/3] Matching amplitude (peak w'' RMS) ...\n');
result = findAmplitude(StabGrid, piv, in);

%% Report + plot
fprintf('\n================= AmplitudeFinisher =================\n');
fprintf('  matched x/c          : %s\n', num2str(piv.xc(sel)));
fprintf('  spanwise modes N     : %d\n', result.N);
fprintf('  matched A0_fund      : %.6g\n', result.A0);
fprintf('  (best swept A0       : %.6g,  J = %.4g)\n', result.A0sweptBest, result.J);
fprintf('=====================================================\n');

save(fullfile(here, 'AmplitudeFinisher_match.mat'), 'result', 'piv', 'sel', 'in');
plot_match(result, piv, sel, here);


% ---------------------------------------------------------------------
function plot_match(result, piv, sel, outdir)
% Overlay HNS vs PIV w'-RMS profiles at the matched stations. Uses a grid
% layout so it stays readable whether the window is 2 stations or 12
% (multi-station collapse test).
hns = result.hns; nS = numel(sel);
nc = min(nS, 5); nr = ceil(nS/nc);
fig = figure('Position', [40 40 300*nc 340*nr]);
for k = 1:nS
    i = sel(k); subplot(nr, nc, k);
    [~, ix] = min(abs(hns.x - piv.x(i)));
    hn = interp1(hns.y, hns.rmsFull(:, ix), piv.y{i}, 'linear', 'extrap');
    if isfield(piv,'rmsFull') && ~isempty(piv.rmsFull{i}); g = piv.rmsFull{i}; else; g = piv.rms{i}; end
    plot(g, piv.y{i}*1e3, 'ko-', 'MarkerSize', 3); hold on;
    plot(hn, piv.y{i}*1e3, 'r-', 'LineWidth', 2);
    if mod(k-1,nc)==0; ylabel('y [mm]'); end
    if k > nS-nc; xlabel('w'' RMS [m/s]'); end
    title(sprintf('x/c = %g%%', piv.xc(i))); grid on; ylim([0 4]);
    if k==1; legend('PIV', 'HNS', 'Location', 'northeast', 'FontSize', 7); end
end
sgtitle(sprintf('AmplitudeFinisher: A0 = %.4g (N=%d)', result.A0, result.N));
saveas(fig, fullfile(outdir, 'AmplitudeFinisher_match.png'));
fprintf('  saved: %s\n', fullfile(outdir, 'AmplitudeFinisher_match.png'));
end
