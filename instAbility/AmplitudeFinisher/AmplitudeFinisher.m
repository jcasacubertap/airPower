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

%% 1. Base flow -> StabGrid  (problem-specific: flat | curved)
fprintf('\n[1/4] Building StabGrid from base flow ...\n');
[StabGrid, in] = buildStabGrid(in);

%% 2. PIV target (experimental .mat placed in AmplitudeFinisher/io/)
fprintf('[2/4] Loading PIV: %s ...\n', in.pivName);
piv = loadPIVwPrime(in.pivFile, in);
in.validation.nModes = piv.nModes;   % HNS diagnostics match the PIV modes present

%% 3. Linear check (ALWAYS) -> figure + implied A0_fund that seeds the search
% One LINEAR HNS, its w'-RMS anchored to PIV at a single station
% (in.linearCheck.xcAnchor). A linear solution scales uniformly with the inlet
% amplitude, so the one anchor scale applies everywhere: the up/downstream
% comparison shows whether the LINEAR growth rate matches the experiment, and
% the implied inlet A0_fund is the physical linear estimate — which then seeds
% the nonlinear sweep below. Runs on a domain covering ALL PIV stations.
xEndL = max(piv.x) / in.match.bufferFrac;
[SGlin, tL] = truncateStabGrid(StabGrid, xEndL, in);
fprintf('[3/4] Linear check: one linear HNS anchored at x/c = %g%% (domain nx %d -> %d) ...\n', ...
        in.linearCheck.xcAnchor, tL.nxOld, tL.nxNew);
lc = linear_check(SGlin, piv, in);
fprintf('\n=============== linear check ===============\n');
fprintf('  anchor x/c            : %g%%\n', in.linearCheck.xcAnchor);
fprintf('  implied A0_fund (lin) : %.6g   <- seeds the nonlinear sweep\n', lc.A0implied);
fprintf('============================================\n');
plot_linear_check(lc, piv, in, here);
fprintf('  saved: %s\n', fullfile(here, 'AmplitudeFinisher_linearCheck.png'));

%% 4. Nonlinear amplitude match over the chordwise window (seeded by lc.A0implied)
sel = find(piv.xc >= in.match.xcWindow(1) & piv.xc <= in.match.xcWindow(2));
if isempty(sel)
    error('AmplitudeFinisher:emptyWindow', ...
          'No PIV stations with x/c in [%g %g]. Available x/c: %s', ...
          in.match.xcWindow(1), in.match.xcWindow(2), num2str(piv.xc));
end
in.match.stationIdx = sel;
xEnd = piv.x(sel(end)) / in.match.bufferFrac;
[SGmatch, tinfo] = truncateStabGrid(StabGrid, xEnd, in);
fprintf('[4/4] Matching amplitude (peak w'' RMS) ; x/c = %s ; domain nx %d -> %d\n', ...
        num2str(piv.xc(sel)), tinfo.nxOld, tinfo.nxNew);
result = findAmplitude(SGmatch, piv, in, lc.A0implied);
StabGrid = SGmatch;   % the (window-truncated) grid actually used by the match

%% Report + plot
fprintf('\n================= AmplitudeFinisher =================\n');
fprintf('  matched x/c          : %s\n', num2str(piv.xc(sel)));
fprintf('  spanwise modes N     : %d\n', result.N);
fprintf('  matched A0_fund      : %.6g\n', result.A0);
fprintf('  match residual J     : %.4g   (0 = perfect fit)\n', result.J);
fprintf('=====================================================\n');

% --- Save all relevant information and fields (metadata + available fields) ---
% Bundles the match outcome (result: A0, sweep, best-swept w'), the full input
% config (in), the PIV target (piv) and matched stations (sel), the base-flow
% grid actually used (StabGrid: base flow + lref/Uref, truncated to the window),
% and a headline `meta` summary. NOTE: there is no full StabRes here — result.A0
% is the interpolated ratio->1 crossing and no solve is run at it; result.hns is
% the extracted w' at A0sweptBest. To post-process the matched solution, run the
% matching DeHNSSo caller (in.caller: sweptwing_flat.m/m3j.m) with A0_fund=result.A0 on the full
% domain (verbatim, no /2).
meta = struct( ...
    'A0_fund',       result.A0, ...          % <-- plug into DeHNSSo caller VERBATIM (no /2)
    'J',             result.J, ...           % match residual (0 = perfect fit over the window)
    'N',             result.N, ...
    'xcWindow',      in.match.xcWindow, ...
    'xcMatched',     piv.xc(sel), ...
    'lambda_z',      in.Stab.lambda_z, ...
    'rmsFactor',     in.match.rmsFactor, ...  % convention: sqrt(2) (mode |a|=A/2 -> z-RMS)
    'xOffset',       in.match.xOffset, ...
    'Uref',          StabGrid.Uref, ...
    'lref',          StabGrid.lref, ...
    'pivFile',       in.pivName, ...
    'savedOn',       char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')), ...
    'problem',       in.problem, ...
    'caller',        in.caller, ...
    'note',          sprintf('A0_fund = result.A0 -> DeHNSSo caller %s.m verbatim (no /2); rmsFactor=sqrt(2)', in.caller));

matFile = fullfile(here, 'AmplitudeFinisher_match.mat');
save(matFile, 'result', 'in', 'piv', 'sel', 'StabGrid', 'meta', '-v7.3');
fprintf('  saved: %s\n', matFile);
fprintf('  -> A0_fund = %.6g  (caller line 28, verbatim)\n', meta.A0_fund);

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


% ---------------------------------------------------------------------
function plot_linear_check(lc, piv, in, outdir)
% Overlay PIV vs the amplitude-anchored LINEAR HNS w'-RMS at ALL stations, so
% both up- and down-stream of the anchor are visible. The linear curve is scaled
% by the single anchor factor lc.c (peak-matched to PIV at xcAnchor); if the
% linear growth rate is right it tracks PIV everywhere, otherwise it diverges
% up/downstream — the linear-vs-nonlinear verdict.
hns = lc.hns;  nS = numel(piv.x);
nc = min(nS, 5);  nr = ceil(nS/nc);
fig = figure('Position', [40 40 300*nc 340*nr]);
for k = 1:nS
    subplot(nr, nc, k);
    [~, ix] = min(abs(hns.x - piv.x(k)));
    hn = lc.c * interp1(hns.y, hns.rmsFull(:, ix), piv.y{k}, 'linear', 'extrap');
    if isfield(piv,'rmsFull') && ~isempty(piv.rmsFull{k}); g = piv.rmsFull{k}; else; g = piv.rms{k}; end
    plot(g, piv.y{k}*1e3, 'ko-', 'MarkerSize', 3); hold on;
    plot(hn, piv.y{k}*1e3, 'r-', 'LineWidth', 2);
    if mod(k-1,nc)==0; ylabel('y [mm]'); end
    if k > nS-nc; xlabel('w'' RMS [m/s]'); end
    ttl = sprintf('x/c = %g%%', piv.xc(k));
    if k == lc.ianchor; ttl = [ttl ' (anchor)']; end
    title(ttl); grid on; ylim([0 4]);
    if k==1; legend('PIV', 'linear HNS', 'Location', 'northeast', 'FontSize', 7); end
end
sgtitle(sprintf('Linear check: anchored at x/c = %g%%,  implied A0_{fund} = %.4g  (N=%d)', ...
        in.linearCheck.xcAnchor, lc.A0implied, lc.N));
saveas(fig, fullfile(outdir, 'AmplitudeFinisher_linearCheck.png'));
fprintf('  saved: %s\n', fullfile(outdir, 'AmplitudeFinisher_linearCheck.png'));
end
