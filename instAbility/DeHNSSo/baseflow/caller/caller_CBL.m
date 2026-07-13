
clc; clear all; format long
caller_dir = fileparts(mfilename('fullpath'));
rootdir    = fullfile(caller_dir, '..', '..');     % refactor repo root
addpath(genpath(rootdir));

%% ===== USER INPUTS =====

% Case: Zanus/Bertolotti 2017 PSE CPG benchmark, Ma=0.5, Te=206 K, adiabatic flat plate.
% Mirror of tests/configs/case_M500_F40.m.

% Flow parameters
BL.Ma       = 0.50;
BL.Pr       = 0.72;
BL.gamma    = 1.4;

% Boundary layer type
BL.beta_H   = 0;            % Hartree parameter (0 = flat plate)
BL.We_ratio = 0;            % We/Ue crossflow ratio (0 = 2D)
BL.wall_bc  = 'adiabatic';
BL.Tw       = 1.0;          % wall T/T_inf (isothermal only)

% Similarity solver settings
BL.eta_inf  = 100;
BL.npts     = 501;

% Freestream thermodynamic state
T_e      = 206;             % [K]
p_e      = 1500;            % [Pa]
R_gas    = 287.1;           % [J/(kg K)] air
mu_ref   = 1.716e-5;        % Sutherland reference μ
T_ref_mu = 273.15;          % Sutherland reference T
S_suth   = 110.4;           % [K]

% Non-dim domain (Blasius lref convention: R_start = Re_0).  Zanus/Bertolotti M05 data
% spans R ~ 643 to 1381 — start upstream, end past branch II + buffer.
R_start  = 500;
R_end    = 1700;
Hn       = 300;             % wall-normal extent in lref_0 units
yi_n     = 10;              % Malik mapping median in lref_0 units

BL.nx    = 800;
BL.ny    = 81;

%% ===== DERIVED QUANTITIES =====

BL.Sstar = S_suth / T_e;
a_e      = sqrt(BL.gamma * R_gas * T_e);
V_inf    = BL.Ma * a_e;
mu_e     = mu_ref * (T_e/T_ref_mu)^1.5 * (T_ref_mu + S_suth) / (T_e + S_suth);
rho_e    = p_e / (R_gas * T_e);
BL.nu    = mu_e / rho_e;

lref_0   = R_start * BL.nu / V_inf;
BL.S     = R_start^2 * BL.nu / V_inf;
BL.L     = R_end^2   * BL.nu / V_inf - BL.S;
BL.H     = Hn   * lref_0;
BL.y_i   = yi_n * lref_0;

BL.Ue    = ones(1, BL.nx) * V_inf;
BL.We    = 0;

Uref = V_inf;
lref = sqrt(BL.S * BL.nu / Uref);
Re_0 = Uref * lref / BL.nu;

%% Similarity solver

[sol, BL] = compressibleFSC(BL);

%% Build physical BL fields on the Chebyshev grid, then non-dimensionalize

nx = BL.nx;  ny = BL.ny;
x_phys = linspace(BL.S, BL.S + BL.L, nx);

% Malik-Chebyshev grid in physical y [m] — only y_phys is needed downstream
[yvec, DM] = chebdif(ny, 2);
[y_phys, ~, ~] = MappingMalik(BL.H, BL.y_i, yvec', DM(:,:,1), DM(:,:,2));

% Similarity profiles (index 1 = wall, in physical eta)
eta_sim = BL.eta;

% Build dimensional fields by interpolating similarity profiles at each station
u   = zeros(ny, nx);  v   = zeros(ny, nx);  w   = zeros(ny, nx);
rho = zeros(ny, nx);  T   = zeros(ny, nx);

for i = 1:nx
    xi   = x_phys(i);
    Ue_i = BL.Ue(i);
    % Similarity variable consistent with f''' + f*f'' = 0 convention
    eta_local = y_phys * sqrt(Ue_i / (BL.nu * xi));

    U_1d      = interp1(eta_sim, BL.fp,       eta_local, 'spline', BL.fp(end));
    W_1d      = interp1(eta_sim, BL.g,        eta_local, 'spline', BL.g(end));
    f_1d      = interp1(eta_sim, BL.f,        eta_local, 'spline', BL.f(end));
    rho_1d    = interp1(eta_sim, BL.rho_prof, eta_local, 'spline', BL.rho_prof(end));
    T_1d      = interp1(eta_sim, BL.theta,    eta_local, 'spline', BL.theta(end));
    etaLD_1d  = interp1(eta_sim, BL.eta_LD,   eta_local, 'spline', BL.eta_LD(end));

    u(:,i)   = Ue_i * U_1d;
    w(:,i)   = BL.We * W_1d;
    % V from Lees-Dorodnitsyn back-transformation (general Ma):
    %   v = sqrt(nu*Ue/(2*x)) * (eta_LD*f' - f) / rho
    % where eta_LD and f are both in the LD coordinate system.
    v(:,i)   = sqrt(BL.nu * Ue_i / (2 * xi)) .* (etaLD_1d .* U_1d - f_1d) ./ rho_1d;
    rho(:,i) = rho_1d;
    T(:,i)   = T_1d;
end

%% Assemble BF struct (non-dim physical fields only — gridgen builds D1/D2 + derivatives)

% Reference scales (Blasius lref convention)
BF.Ma    = BL.Ma;
BF.Pr    = BL.Pr;
BF.gamma = BL.gamma;
BF.Sstar = BL.Sstar;
BF.Uref  = Uref;
BF.lref  = lref;
BF.nu    = BL.nu;
BF.Re    = Re_0;

% Non-dim grid
BF.X   = x_phys / lref;       % [1 x nx]
BF.Y   = (y_phys / lref).';   % [ny x 1]
BF.H   = BL.H   / lref;
BF.y_i = BL.y_i / lref;

% Non-dim fields
BF.U   = u   / Uref;
BF.V   = v   / Uref;
BF.W   = w   / Uref;
BF.rho = rho;                 % already rho/rho_inf
BF.T   = T;                   % already T/T_inf

%% Save

outdir = fullfile(caller_dir, '..', 'output');
if ~exist(outdir, 'dir'); mkdir(outdir); end
save(fullfile(outdir, 'BF_compressible.mat'), 'BF');

fprintf('\nSaved:\n');
fprintf('  BF -> baseflow/output/BF_compressible.mat\n');
fprintf('  nx=%d, ny=%d, Re_0=%.2f, Ma=%.3f, beta_H=%.3f, We/Ue=%.2f\n', ...
    nx, ny, Re_0, BL.Ma, BL.beta_H, BL.We_ratio);

%% Postprocessing (BF profiles at inflow / mid station)

figure('Position', [100 100 1400 350]);
subplot(1,5,1); plot(BF.U(:,1),             BF.Y); xlabel('U/U_e');         ylabel('y/l_{ref}'); title('U');
subplot(1,5,2); plot(BF.W(:,1),             BF.Y); xlabel('W/W_e');         title('W');
subplot(1,5,3); plot(BF.T(:,1),             BF.Y); xlabel('T/T_\infty');    title('T');
subplot(1,5,4); plot(BF.rho(:,1),           BF.Y); xlabel('\rho/\rho_\infty'); title('\rho');
subplot(1,5,5); plot(BF.V(:,round(nx/2)),   BF.Y); xlabel('V');             title('V');
