function R = destreak_feasibility(cfg)
%DESTREAK_FEASIBILITY  Frozen-base-flow inner one-shot for passive de-streaking.
%
%   R = destreak_feasibility(cfg)
%
%   Builds the wall-relief -> induced-roller operator V (and streak-change S) by
%   finite-difference geometry sensitivity through gridgen, forms the sealed-
%   Poisson ceiling P, and answers the decisive feasibility questions (brief S6):
%     (a) singular values of V         -> is the reachable roller subspace rich?
%     (b) overlap of range(V) with the productive drain direction d = -U'*uhat
%                                       -> can a spanwise-invariant relief drain at all?
%     (c) one least-squares solve for c*, achieved transfer T=-I2 over the window,
%         and efficiency eta vs the sealed-Poisson ceiling.
%
%   ONE frozen base flow, ONE inner solve per relief mode. No outer loop, no
%   base-flow regeneration (deferred OpenFOAM coupling). Uses only the DeHNSSo
%   tree inside OptimalDecayNumerics.

if nargin < 1; cfg = struct(); end
def = struct('n_eta',50,'n_xi',800,'lu_mode','none','K',6);
fn = fieldnames(def); for i=1:numel(fn); if ~isfield(cfg,fn{i}); cfg.(fn{i})=def.(fn{i}); end; end

here = fileparts(mfilename('fullpath'));  addpath(here);
addpath(genpath(fullfile(here,'..','DeHNSSo','src')));
addpath(genpath(fullfile(here,'..','DeHNSSo','gridgen','src')));

%% 1. Flat frozen base flow
gridfile = fullfile(here, sprintf('flat_grid_%dx%d.mat', cfg.n_eta, cfg.n_xi));
if isfile(gridfile)
    S = load(gridfile);  G = S.StabGrid;
    fprintf('[feas] loaded %s\n', gridfile);
else
    G = destreak_make_flatgrid(cfg.n_eta, cfg.n_xi, gridfile);
end
beta0 = 2*pi * G.lref / 7.5e-3;
fprintf('[feas] ny=%d nx=%d beta0=%.4f\n', G.ny, G.nx, beta0);

%% 1b. Lift the flat grid onto the curvilinear representation (2-D X,Y).
%  DY=0 must reproduce the flat operator: identity metrics, same dyU.
Gc = destreak_perturb_grid(G, zeros(G.ny, G.nx));
rel_dyU = max(abs(Gc.dyU(:) - G.dyU(:))) / max(abs(G.dyU(:)));
fprintf('[feas] curvilinear-lift reproduces flat: rel dyU diff = %.2e\n', rel_dyU);
assert(rel_dyU < 1e-6, 'curvilinear lift does not reproduce flat grid (rel=%.2e)', rel_dyU);

%% 2. Flat streak (homogeneous solution)  u_hom, v_hom  (on the curvilinear base)
ropt.lu_mode = cfg.lu_mode;
[u0, v0] = destreak_run_case(Gc, beta0, ropt);
fprintf('[feas] flat streak: max|u|=%.3e max|v_hom|=%.3e\n', max(abs(u0(:))), max(abs(v0(:))));

%% 3. Relief basis and window
rcfg = struct('K',cfg.K);  if isfield(cfg,'relief'); rcfg = cfg.relief; end
B = destreak_relief_basis(Gc, rcfg);
stn = find(B.win);  ns = numel(stn);  K = size(B.gmat,2);
fprintf('[feas] K=%d relief modes, window = %d stations\n', K, ns);

%% 4. Geometry-sensitivity sweep: columns of V (roller) and S (streak change)
ny = Gc.ny;
Vfull = zeros(ny, Gc.nx, K);
Sfull = zeros(ny, Gc.nx, K);
for k = 1:K
    Gk = destreak_perturb_grid(Gc, B.DY{k});
    [uk, vk] = destreak_run_case(Gk, beta0, ropt);
    % phase-align to the flat run (ILST fixes an arbitrary global phase per run)
    a = sum(conj(uk(:)).*u0(:));  ph = a/abs(a);
    uk = uk*ph;  vk = vk*ph;
    Sfull(:,:,k) = (uk - u0) / B.eps;
    Vfull(:,:,k) = (vk - v0) / B.eps;
    fprintf('  [feas] mode %d/%d: max|dV|/eps=%.3e max|dS|/eps=%.3e\n', ...
            k, K, max(abs(Vfull(:,:,k)),[],'all'), max(abs(Sfull(:,:,k)),[],'all'));
end

%% 5. Sealed-Poisson ceiling operator P over the window
[applyP, Vopt_hom] = destreak_Popt(Gc, beta0, stn, u0);   % P u_hom  [ny x ns]

%% 6. Assemble the least-squares system   A c ~= r  over (y, station) in window
% residual = v_hat(c) - P u_hat(c) = (v0 - P u0) + (V - P S) c
% => A = V - P S,  r = P u0 - v0
flat = @(M) reshape(M, [], 1);
r  = flat(Vopt_hom - v0(:,stn));                 % [ny*ns x 1]
A  = zeros(ny*ns, K);
for k = 1:K
    PSk = applyP(Sfull(:,:,k));                   % [ny x ns]
    A(:,k) = flat(Vfull(:,stn,k) - PSk);
end

%% 6b. Metric weights (diagonal). Energy: trapz(y)*dx. g-weighted: * |drain|.
x = Gc.x(:).';  dxv = gradient(x);
wq = zeros(ny, ns);
for idx = 1:ns
    i = stn(idx);
    wy = trap_weights(Gc.Y(:,i));
    wq(:,idx) = abs(wy) * abs(dxv(i));
end
wqv   = flat(wq);
dvec  = flat(-Gc.dyU(:,stn).*u0(:,stn));          % productive drain roller direction
gw    = abs(dvec)/max(abs(dvec));
W_E   = wqv;                                       % energy metric
W_g   = wqv .* gw;                                 % g-weighted (productive) metric

%  The relief is REAL: h(x)=sum c_k g_k(x), c_k in R. The induced roller is a
%  complex field linear in real c, so the least-squares is a REAL problem over
%  the stacked (Re;Im) of the complex operator.
stack = @(M) [real(M); imag(M)];

%% 7. Real-c least-squares solve for c* (both metrics) and reachability
solve_c = @(W) (sqrt([W;W]).*stack(A)) \ (sqrt([W;W]).*stack(r));
c_E = solve_c(W_E);
c_g = solve_c(W_g);

% Weighted R^2 of the roller-correction fit = fraction of the required roller
% change r = (P u_hom - v_hom) that a REAL spanwise-invariant relief can produce.
r2 = @(c,W) 1 - sum(([W;W]).*(stack(A*c - r)).^2) / sum(([W;W]).*stack(r).^2);
R2_E = r2(c_E, W_E);
R2_g = r2(c_g, W_g);

% Richness of the real reachable roller subspace (singular values of V under
% real combinations): stack Re/Im of the weighted V columns.
sWs    = sqrt([W_E; W_E]);
Vstack = sWs .* stack(reshape(Vfull(:,stn,:), ny*ns, K));
sv     = svd(Vstack);

%% 8. Achieved fields and transfer T=-I2 over the window (draining sign: T<0)
uA = u0;  vA = v0;
for k=1:K; uA = uA + c_g(k)*Sfull(:,:,k); vA = vA + c_g(k)*Vfull(:,:,k); end
T_flat = transfer(Gc, u0, v0,          stn);            % streak's own (sustaining) roller
T_ceil = transfer(Gc, u0, embed(Vopt_hom,stn,Gc), stn);% sealed-Poisson ceiling (draining)
T_ach  = transfer(Gc, uA, vA,          stn);            % achieved design
% Progress from the flat (sustaining) state toward the draining ceiling.
progress = (sum(T_flat) - sum(T_ach)) / (sum(T_flat) - sum(T_ceil));

%% 9. Report + save
fprintf('\n==================== FEASIBILITY ====================\n');
fprintf(' sing. values of V (real reachable): %s\n', mat2str(sv(:)',3));
fprintf(' roller-correction reach  R^2 (E)  : %.3f   (frac of needed roller reachable)\n', R2_E);
fprintf(' roller-correction reach  R^2 (g)  : %.3f\n', R2_g);
fprintf(' c* (g-weighted, REAL)             : %s\n', mat2str(c_g(:)',3));
fprintf(' sum T_flat  (flat self-roller)    : %+.4e   [+ = sustains streak]\n', sum(T_flat));
fprintf(' sum T_ceil  (sealed-Poisson opt)  : %+.4e   [- = drains streak]\n', sum(T_ceil));
fprintf(' sum T_ach   (design c*)           : %+.4e\n', sum(T_ach));
fprintf(' progress flat->ceiling            : %.3f   (1 = full de-streaking)\n', progress);
fprintf('====================================================\n');

R = struct('G',Gc,'beta0',beta0,'B',B,'stn',stn, ...
           'u0',u0,'v0',v0,'Vfull',Vfull,'Sfull',Sfull,'Vopt_hom',Vopt_hom, ...
           'A',A,'r',r,'W_E',W_E,'W_g',W_g,'sv',sv, ...
           'R2_E',R2_E,'R2_g',R2_g,'c_E',c_E,'c_g',c_g, ...
           'T_flat',T_flat,'T_ceil',T_ceil,'T_ach',T_ach,'progress',progress);
save(fullfile(here, sprintf('feas_%dx%d.mat', cfg.n_eta, cfg.n_xi)), 'R', '-v7.3');
end

%% ---- helpers ----
function w = trap_weights(Y)
% trapezoidal quadrature weights for arbitrary (monotone) nodes Y
Ys = flipud(Y(:));  d = diff(Ys);  n = numel(Ys);
w = zeros(n,1);  w(1)=d(1)/2;  w(end)=d(end)/2;  w(2:end-1)=(d(1:end-1)+d(2:end))/2;
w = flipud(w);
end

function T = transfer(G, uhat, vhat, stn)
% T(x) = -I2 = -int U' Re(conj(uhat) vhat) dy.
% Sign convention (verified against the sealed-Poisson optimal):
%   T < 0  drains the streak (de-streaking, = the ceiling's sign)
%   T > 0  sustains/amplifies the streak (normal lift-up)
ns = numel(stn);  T = zeros(1,ns);
for idx=1:ns
    i = stn(idx);
    integ = G.dyU(:,i) .* real(conj(uhat(:,i)).*vhat(:,i));
    w = trap_weights(G.Y(:,i));
    T(idx) = -sum(w .* integ);
end
end

function Vfullcol = embed(Vwin, stn, G)
% place a [ny x ns] window field into a full [ny x nx] array (zeros elsewhere)
Vfullcol = zeros(G.ny, G.nx);
Vfullcol(:,stn) = Vwin;
end
