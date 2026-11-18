function B = destreak_relief_basis(G, cfg)
%DESTREAK_RELIEF_BASIS  Spanwise-invariant wall-relief basis g_k(x) and taper phi(eta).
%
%   B = destreak_relief_basis(G, cfg)
%
%   Builds a small basis of localized bumps g_k(x) (Gaussians) spanning the
%   design window [xa,xb], plus a wall-normal taper phi(eta) (=1 at the wall,
%   ->0 at the freestream). The relief is h(x) = sum_k c_k g_k(x); the k-th mesh
%   displacement is  DY_k(eta,x) = eps * phi(eta) * g_k(x).
%
%   cfg fields (all optional):
%     .K          number of bumps                     (default 6)
%     .win_frac   [lo hi] window as arclength frac     (default [0.30 0.72])
%     .eps_frac   bump height / domain height          (default 5e-3)
%     .s0         taper decay (frac of BL height)      (default 0.18)
%     .sig_frac   bump width as frac of window span    (default 1/(2.5*K))
%
%   Returns B with:
%     .gmat   [nx x K]   g_k(x), each max-normalized to 1
%     .phi    [ny x 1]   wall taper phi(eta)
%     .eps    scalar     physical displacement amplitude
%     .DY     1xK cell   DY_k = eps*phi*g_k'   [ny x nx]
%     .win    logical 1xnx   window mask (design stations)
%     .xc     1xK        bump centres (in x units)

if nargin < 2; cfg = struct(); end
def = struct('K',6,'win_frac',[0.30 0.72],'eps_frac',5e-4,'s0',0.18);
fn = fieldnames(def);
for i=1:numel(fn); if ~isfield(cfg,fn{i}); cfg.(fn{i})=def.(fn{i}); end; end
K = cfg.K;
if ~isfield(cfg,'sig_frac'); cfg.sig_frac = 1/(2.5*K); end

x   = G.x(:).';               % 1 x nx  (streamwise, = xi + Re)
mid = round(G.nx/2);
Yc  = G.Y(:,mid);             % representative column
Ywall = Yc(G.ny);  Ytop = Yc(1);
Hcol  = Ytop - Ywall;

%% Wall-normal taper phi(eta): 1 at wall, Gaussian decay to freestream
s   = (Yc - Ywall) / Hcol;    % 0 at wall (row ny), 1 at top (row 1)
phi = exp(-(s / cfg.s0).^2);
phi(end) = 1;                 % enforce exactly 1 at the wall node

%% Window and bump centres
xn   = (x - x(1)) / (x(end) - x(1));       % 0..1 arclength
win  = xn >= cfg.win_frac(1) & xn <= cfg.win_frac(2);
span = cfg.win_frac(2) - cfg.win_frac(1);
pad  = 0.10 * span;
xc_n = linspace(cfg.win_frac(1)+pad, cfg.win_frac(2)-pad, K);
sig_n = cfg.sig_frac * span;
% map normalized centres/width back to x units
xc  = x(1) + xc_n * (x(end)-x(1));
sig = sig_n * (x(end)-x(1));

%% Bumps g_k(x) — Gaussians, max-normalized, ~0 at inflow (protects ILST)
gmat = zeros(G.nx, K);
for k = 1:K
    g = exp(-((x - xc(k))/sig).^2).';
    g = g / max(g);
    gmat(:,k) = g;
end

eps_amp = cfg.eps_frac * Hcol;
DY = cell(1,K);
for k = 1:K
    DY{k} = eps_amp * (phi * gmat(:,k).');   % [ny x nx]
end

% Guard: displacement at the inflow station must be negligible (frozen ILST).
inflow_disp = max(cellfun(@(d) max(abs(d(:,1))), DY));
if inflow_disp > 1e-6 * eps_amp
    warning('destreak:inflowDisp', ...
        'relief has %.2e displacement at inflow — ILST inflow may drift', inflow_disp);
end

B = struct('gmat',gmat,'phi',phi,'eps',eps_amp,'DY',{DY}, ...
           'win',win,'xc',xc,'sig',sig,'cfg',cfg);
end
