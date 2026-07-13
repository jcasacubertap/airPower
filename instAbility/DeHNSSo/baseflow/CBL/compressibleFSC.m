function [sol, BL] = compressibleFSC(BL)
%COMPRESSIBLEFSC  Compressible Falkner-Skan-Cooke boundary layer solver.
%
%   [sol, BL] = compressibleFSC(BL)
%
%   Solves the compressible similarity equations in Lees-Dorodnitsyn variables
%   for a swept wedge (Falkner-Skan-Cooke with energy equation):
%
%     (C f'')' + f f'' + beta_H (T - f'^2) = 0      (streamwise momentum)
%     (C g')'  + f g'                       = 0      (crossflow momentum)
%     (C/Pr theta')' + f theta' + Ec C (f''^2 + g'^2) = 0   (energy)
%
%   where  C = rho*mu = mu(T)/T   (Chapman-Rubesin parameter),
%          T = theta   (non-dim temperature T/T_inf),
%          rho = 1/T   (ideal gas, constant pressure across BL),
%          Ec = (gamma-1)*Ma^2   (Eckert number).
%
%   Sutherland viscosity:  mu(T) = T^(3/2) (1+S*)/(T+S*),  S* = S/T_inf.
%
%   BCs:  f(0)=0, f'(0)=0, f'(inf)=1    (no-slip, edge velocity)
%         g(0)=0, g(inf)=1               (no-slip, edge crossflow)
%         theta(inf)=1                    (edge temperature)
%         theta(0)=Tw/T_inf (isothermal) or theta'(0)=0 (adiabatic)
%
%   For beta_H=0, g=0: recovers compressible Blasius.
%   For Ma->0 (C->1, T->1): recovers incompressible FSC.
%
%   Input struct BL fields:
%     Ma       — freestream Mach number
%     Pr       — Prandtl number (default 0.72)
%     gamma    — ratio of specific heats (default 1.4)
%     Sstar    — Sutherland constant S/T_inf (default 110.4/288)
%     beta_H   — Hartree parameter (0 = flat plate, 1 = stagnation)
%     wall_bc  — 'adiabatic' or 'isothermal' (default 'adiabatic')
%     Tw       — wall temperature T_wall/T_inf (required if isothermal)
%     eta_inf  — outer boundary for similarity variable (default 20)
%     npts     — number of BVP collocation points (default 401)
%
%   Output:
%     sol      — bvp5c solution structure
%     BL       — input struct augmented with solution profiles:
%       .eta, .f, .fp, .fpp, .g, .gp, .theta, .thetap
%       .rho, .mu, .muT, .V   (derived quantities)

%% Defaults
if ~isfield(BL,'Pr');      BL.Pr = 0.72; end
if ~isfield(BL,'gamma');   BL.gamma = 1.4; end
if ~isfield(BL,'Sstar');   BL.Sstar = 110.4/288; end
if ~isfield(BL,'beta_H');  BL.beta_H = 0; end
if ~isfield(BL,'wall_bc'); BL.wall_bc = 'adiabatic'; end
if ~isfield(BL,'eta_inf'); BL.eta_inf = 20; end
if ~isfield(BL,'npts');    BL.npts = 401; end
if ~isfield(BL,'Tw');      BL.Tw = 1; end
if ~isfield(BL,'We_ratio'); BL.We_ratio = 0; end

Ma   = BL.Ma;
Pr   = BL.Pr;
ga   = BL.gamma;
Ss   = BL.Sstar;
bH   = BL.beta_H;
Ec   = (ga - 1) * Ma^2;
Wr2  = BL.We_ratio^2;        % crossflow dissipation scaling
is_adiabatic = strcmpi(BL.wall_bc, 'adiabatic');

fprintf('Compressible FSC solver: Ma=%.3f, beta_H=%.3f, %s wall\n', ...
    Ma, bH, BL.wall_bc);

%% State vector:  y = [f, f', f'', g, g', theta, theta']
%  Indices:        1   2    3    4   5     6       7

    function C = chapman(T)
        % C = rho*mu = mu(T)/T,  with mu = T^(3/2)(1+S*)/(T+S*)
        C = T.^(1/2) .* (1+Ss) ./ (T+Ss);
    end

    function Cp = dchapman(T)
        % dC/dT
        mu  = T.^(3/2) .* (1+Ss) ./ (T+Ss);
        muT = mu./T .* (3/2 - T./(T+Ss));
        Cp  = (muT.*T - mu) ./ T.^2;  % d(mu/T)/dT = (muT*T - mu)/T^2
    end

    function dydt = odefun(~, y)
        fp  = y(2);  fpp = y(3);
        gp  = y(5);
        th  = y(6);  thp = y(7);

        C  = chapman(th);
        Cp = dchapman(th);

        % (C f'')' = C f''' + C' theta' f''
        % => f''' = [-y(1)*f'' - beta_H*(th - fp^2) - Cp*thp*fpp] / C
        fppp = (-y(1)*fpp - bH*(th - fp^2) - Cp*thp*fpp) / C;

        % (C g')' = C g'' + C' theta' g'
        % => g'' = [-y(1)*g' - Cp*thp*gp] / C
        gpp = (-y(1)*gp - Cp*thp*gp) / C;

        % (C/Pr theta')' = C/Pr theta'' + Cp/Pr thp^2
        % Dissipation: Ec C (f''^2 + We^2 g'^2),  scaled by We_ratio^2
        % => theta'' = [-Pr*y(1)*thp - Pr*Ec*C*(fpp^2 + Wr2*gp^2) - Cp*thp^2] / C
        thpp = (-Pr*y(1)*thp - Pr*Ec*C*(fpp^2 + Wr2*gp^2) - Cp*thp^2) / C;

        dydt = [fp; fpp; fppp; gp; gpp; thp; thpp];
    end

    function res = bcfun(ya, yb)
        res = zeros(7, 1);
        res(1) = ya(1);        % f(0) = 0
        res(2) = ya(2);        % f'(0) = 0
        res(3) = yb(2) - 1;    % f'(inf) = 1
        res(4) = ya(4);        % g(0) = 0
        res(5) = yb(4) - 1;    % g(inf) = 1  (if no crossflow, g=0 everywhere → set below)
        res(6) = yb(6) - 1;    % theta(inf) = 1
        if is_adiabatic
            res(7) = ya(7);    % theta'(0) = 0  (adiabatic)
        else
            res(7) = ya(6) - BL.Tw;  % theta(0) = Tw  (isothermal)
        end
    end

%% Initial guess
eta_grid = linspace(0, BL.eta_inf, BL.npts);

% Adiabatic recovery temperature as initial guess for theta
Tr = 1 + sqrt(Pr) * Ec / 2;

yinit = @(t) [t - log(t+1);           % f
              1 - 1/(t+1);             % f'
              1/(t+1)^2;               % f''
              1 - exp(-t);             % g
              exp(-t);                 % g'
              1 + (Tr-1)*exp(-t);      % theta (decays from Tr at wall to 1)
              -(Tr-1)*exp(-t)];        % theta'

solinit = bvpinit(eta_grid, yinit);

%% Solve
options = bvpset('AbsTol', 1e-8, 'RelTol', 1e-6, 'NMax', 10000);
warning('off', 'MATLAB:bvp5c:RelTolNotMet');
sol = bvp5c(@odefun, @bcfun, solinit, options);
warning('on', 'MATLAB:bvp5c:RelTolNotMet');

%% Extract profiles on a fine grid
eta_LD = linspace(0, BL.eta_inf, BL.npts)';
Y      = deval(sol, eta_LD)';

BL.eta_LD = eta_LD;          % Lees-Dorodnitsyn coordinate
BL.f      = Y(:,1);
BL.fp     = Y(:,2);          % U/Ue
BL.fpp    = Y(:,3);
BL.g      = Y(:,4);          % W/We
BL.gp     = Y(:,5);
BL.theta  = Y(:,6);          % T/T_inf
BL.thetap = Y(:,7);

% Derived quantities
BL.rho_prof = 1 ./ BL.theta;                                     % ρ/ρ∞
BL.mu_prof  = BL.theta.^(3/2) .* (1+Ss) ./ (BL.theta + Ss);     % μ/μ∞
BL.muT_prof = BL.mu_prof./BL.theta .* (3/2 - BL.theta./(BL.theta+Ss));  % dμ/dT

% Physical similarity variable:  dη_phys/dη_LD = √2 / ρ = √2 T
% Use Simpson's rule (4th-order) for higher accuracy
integrand = sqrt(2) * BL.theta;
n = numel(eta_LD);
eta_phys = zeros(n, 1);
for j = 2:n
    h = eta_LD(j) - eta_LD(j-1);
    % Simpson via bvp5c midpoint interpolation
    mid = (eta_LD(j-1) + eta_LD(j)) / 2;
    Ymid = deval(sol, mid);
    f_mid = sqrt(2) * Ymid(6);  % theta at midpoint
    eta_phys(j) = eta_phys(j-1) + h/6 * (integrand(j-1) + 4*f_mid + integrand(j));
end
BL.eta = eta_phys;            % physical η (use this for StabGrid)

% Wall values
BL.Tw_actual = BL.theta(1);
BL.fpp_wall  = BL.fpp(1);
BL.gp_wall   = BL.gp(1);

fprintf('  Wall: T_w/T_inf = %.4f,  f''''(0) = %.6f,  g''(0) = %.6f\n', ...
    BL.Tw_actual, BL.fpp_wall, BL.gp_wall);

end
