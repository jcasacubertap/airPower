function [StabRes] = CLST(j, i, StabGrid, StabRes, bc_top, Opt)
%CLST  Compressible LST inflow condition for mode j.
%
%   StabRes = CLST(j, i, StabGrid, StabRes, bc_top, Opt)
%
%   Solves the compressible linear stability eigenvalue problem (spatial)
%   at station i for variables (rho, u, v, w, T).  Given omega and beta,
%   finds the streamwise wavenumber alpha and corresponding eigenfunctions.
%
%   Matrix assembly is delegated to CLST_Matrix (one call per grid point).
%   Ideal gas EOS: p = rho*T/(gamma*Ma^2), e = cv*T, Sutherland viscosity.
%
%   Inputs:
%     j         — mode index
%     i         — streamwise station (typically 1 for inflow)
%     StabGrid  — grid structure (U, W, rho, T, D1, D2, dzdy, d2zdy2,
%                 Re, Ma, Pr, gamma, Sstar)
%     StabRes   — stability results (updated with rho, u, v, w, T, phi)
%     bc_top    — freestream BC: 'H_DR' (default) or 'H_NM'
%     Opt       — options struct (Opt.plot = 'on'/'off')

%% Parameters
ny    = StabGrid.ny;
Re    = StabGrid.Re;
Ma    = StabGrid.Ma;
Pr    = StabGrid.Pr;
ga    = StabGrid.gamma;
Sstar = StabGrid.Sstar;
Ec    = (ga - 1) * Ma^2;

omega = StabRes.omegavec(j);
beta  = StabRes.betavec(j);

if ~exist('bc_top','var'); bc_top = 'H_DR'; end
if ~exist('Opt','var');    Opt = struct('plot','on'); end
if ~isfield(Opt,'verbose'); Opt.verbose = 'on'; end

if strcmpi(Opt.verbose,'on')
    fprintf(' Calculating IC for nonzero modes with CLST: ')
end

%% Base flow at station i
U    = StabGrid.U(:,i);
W    = StabGrid.W(:,i);
Rb   = StabGrid.rho(:,i);
Tb   = StabGrid.T(:,i);
dyU  = StabGrid.D1 * U;
dyW  = StabGrid.D1 * W;
dyRb = StabGrid.D1 * Rb;
dyTb = StabGrid.D1 * Tb;

d2yU  = StabGrid.D2 * U;
d2yTb = StabGrid.D2 * Tb;

[mu_b, dmuT, d2muT] = sutherland(Tb, Sstar);
mu_y = dmuT .* dyTb;

%% Chebyshev operators
[~, DM] = chebdif(ny, 2);
Dn  = DM(:,:,1);
Dn2 = Dn * Dn;
I5  = eye(5);
Total_Dn  = kron(Dn,  I5);
Total_DDn = kron(Dn2, I5);

dzdy   = StabGrid.dzdy;
d2zdy2 = StabGrid.d2zdy2;

%% Assemble coefficient matrices
zz = 1i;
N5 = 5 * ny;
AAL1 = zeros(N5);  AAL2 = zeros(N5);  AAL3 = zeros(N5);
BBR1 = zeros(N5);  BBR2 = zeros(N5);

gM2 = ga * Ma^2;
Rg  = 1 / gM2;
cv  = Rg / (ga - 1);  % de/dT = cv (ideal gas, nondimensional)
eT  = cv;

for k = 1:ny
    pk = Rb(k) * Tb(k) / gM2;  % base pressure (nondimensional)

    [Lt, Lx, Ly, Lz, Lq, Vxx, Vxy, Vyy, Vxz, Vyz, Vzz] = CLST_Matrix( ...
        Re, Pr, Ec, ...
        Rb(k), dyRb(k), Tb(k), dyTb(k), d2yTb(k), ...
        U(k), dyU(k), d2yU(k), W(k), dyW(k), ...
        mu_b(k), dmuT(k), d2muT(k), mu_y(k), pk, eT, gM2);

    % Row scaling: multiply energy row (5) by Ec for conditioning
    Lt(5,:) = Lt(5,:) * Ec;   Lx(5,:) = Lx(5,:) * Ec;
    Ly(5,:) = Ly(5,:) * Ec;   Lz(5,:) = Lz(5,:) * Ec;
    Lq(5,:) = Lq(5,:) * Ec;
    Vxx(5,:) = Vxx(5,:) * Ec; Vyy(5,:) = Vyy(5,:) * Ec;
    Vzz(5,:) = Vzz(5,:) * Ec;
    Vxy(5,:) = Vxy(5,:) * Ec; Vxz(5,:) = Vxz(5,:) * Ec;
    Vyz(5,:) = Vyz(5,:) * Ec;

    % Spatial EVP: AA*q = alpha*BB*q  (drop alpha^2 terms)
    A1 = -zz*omega*Lt + zz*beta*Lz + Lq - beta^2*Vzz;
    A2 = Ly + zz*beta*Vyz;
    A3 = Vyy;
    B1 = -zz*Lx + beta*Vxz;
    B2 = -zz*Vxy;

    % Apply Malik mapping
    AA1 = A1;
    AA2 = A2 * dzdy(k) + A3 * d2zdy2(k);
    AA3 = A3 * dzdy(k)^2;
    BB1 = B1;
    BB2 = B2 * dzdy(k);

    icc = (k - 1) * 5;
    AAL1(icc+1:icc+5, icc+1:icc+5) = AA1;
    AAL2(icc+1:icc+5, icc+1:icc+5) = AA2;
    AAL3(icc+1:icc+5, icc+1:icc+5) = AA3;
    BBR1(icc+1:icc+5, icc+1:icc+5) = BB1;
    BBR2(icc+1:icc+5, icc+1:icc+5) = BB2;
end

% Global matrices
AA = AAL1 + AAL2 * Total_Dn + AAL3 * Total_DDn;
BB = BBR1 + BBR2 * Total_Dn;

%% Boundary conditions
% Index 1 = freestream (y = H), index ny = wall (y = 0)
% Ordering within each point: [rho, u, v, w, T]

% Freestream: u = v = w = T = 0  ('H_DR') OR d/dy = 0  ('H_NM')
% Wavenumber BC here is on the wall-normal derivative operator Total_Dn.
use_NM_top = strcmpi(bc_top, 'H_NM');
for m = 2:5
    if use_NM_top
        AA(m,:) = Total_Dn(m,:);   % Neumann: dq^/dy = 0 at top
    else
        AA(m,:) = 0;  AA(m,m) = 1; % Dirichlet (default)
    end
    BB(m,:) = 0;
end

% Wall: u = v = w = T = 0  (last 4 rows of last point)
for m = 0:3
    row = N5 - 3 + m;
    AA(row,:) = 0;  AA(row,row) = 1;
    BB(row,:) = 0;
end

%% Solve eigenvalue problem
if isfield(Opt, 'alpha_shift') && ~isempty(Opt.alpha_shift)
    % Mode continuation: eigs with shift near previous alpha (fast, one EV)
    opts_eigs.tol    = 1e-10;
    opts_eigs.maxit  = 300;
    opts_eigs.v0     = ones(N5, 1);
    [Evec_s, alpha] = eigs(AA, BB, 1, Opt.alpha_shift, opts_eigs);
    alpha = alpha;   % scalar eigenvalue
    index = 1;
    Evec  = Evec_s;
else
    % Full spectrum + EVfilter (used at station 1 for mode identification)
    [Evec, P] = eig(AA, BB);
    Eval = diag(P);

    %% Extract v-component for filtering
    v_evec = Evec(3:5:N5, :);

    kd99 = find(U <= 0.99*max(U), 1, 'first');
    d99_val = StabGrid.y(kd99);

    % Optional Mack-mode selector (Ma>1):
    %   Opt.mode_select = 'auto' (default) | 'slow' | 'fast'
    mode_select = 'auto';
    if isfield(Opt,'mode_select') && ~isempty(Opt.mode_select)
        mode_select = Opt.mode_select;
    end
    [alpha, index] = EVfilter(Eval, v_evec, StabGrid.y, omega, beta, d99_val, Ma, mode_select);
end

StabRes.alpha(j,i) = alpha;
if strcmpi(Opt.verbose,'on')
    fprintf('alpha = %.6f %+.6fi\n', real(alpha), imag(alpha));
end

%% Extract and normalise eigenfunction
%  CLST solves in (rho, u, v, w, T) internally.
%  Convert to (u, v, w, p, T) for the HNS system.
%  Pressure perturbation from EOS: p' = (rho'*T_bar + rho_bar*T') / (gamma*Ma^2)
q       = Evec(1:N5, index);
rho_hat = q(1:5:N5);
u_hat   = q(2:5:N5);
v_hat   = q(3:5:N5);
w_hat   = q(4:5:N5);
T_hat   = q(5:5:N5);
p_hat   = (rho_hat .* Tb + Rb .* T_hat) / gM2;

y_int  = linspace(StabGrid.y(1), StabGrid.y(end), 4000);
u_int  = interp1(StabGrid.y, abs(u_hat), y_int, 'spline');
ampltd = max(abs(u_int));

StabRes.u(j,:,i)   = u_hat.'   / ampltd;
StabRes.v(j,:,i)   = v_hat.'   / ampltd;
StabRes.w(j,:,i)   = w_hat.'   / ampltd;
StabRes.p(j,:,i)   = p_hat.'   / ampltd;
StabRes.rho(j,:,i) = rho_hat.' / ampltd;

% Compressible LHS is density form: slot 4 of phi is rho'.
if isfield(StabRes,'T')
    StabRes.T(j,:,i) = T_hat.' / ampltd;
    StabRes.phi(j,:,i) = [StabRes.u(j,:,i),    StabRes.v(j,:,i), ...
                           StabRes.w(j,:,i),    StabRes.rho(j,:,i), ...
                           StabRes.T(j,:,i)];
else
    StabRes.phi(j,:,i) = [StabRes.u(j,:,i),  StabRes.v(j,:,i), ...
                           StabRes.w(j,:,i),  StabRes.rho(j,:,i)];
end

%% Plotting
if ~strcmpi(Opt.plot, 'on'); return; end

figure(1); clf; set(gcf, 'Position', [200 200 1400 350]);
tiledlayout(1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile
% δ99: y where U first drops below 0.99·U_max (works for ascending OR descending y)
kd99_plot    = find(U >= 0.99*max(U), 1, 'last');
if isempty(kd99_plot); kd99_plot = numel(U); end
delta99_plot = StabGrid.y(kd99_plot);
eta_max_plot = min(2 * delta99_plot, max(StabGrid.y));
plot(U, StabGrid.y, 'b-', 'LineWidth', 1.5); hold on;
plot(Tb, StabGrid.y, 'r-', 'LineWidth', 1.5);
plot(Rb, StabGrid.y, 'g-', 'LineWidth', 1.5);
ylim([0, eta_max_plot]);
xlabel('Base flow', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$\eta$', 'Interpreter', 'latex', 'FontSize', 12)
legend('$U$', '$T$', '$\rho$', 'Location', 'best', 'Interpreter', 'latex')
title('Inflow base flow', 'Interpreter', 'latex', 'FontSize', 13)
grid on; box on;

nexttile
phys = ~isinf(Eval) & ~isnan(Eval) & abs(Eval) < 10;
plot(real(Eval(phys)), imag(Eval(phys)), 'k.', 'MarkerSize', 12); hold on;
plot(real(alpha), imag(alpha), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
xlim([real(alpha)-0.2, real(alpha)+0.2]);
ylim([-0.1, 0.05]);
xlabel('$\alpha_r$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$\alpha_i$', 'Interpreter', 'latex', 'FontSize', 12)
legend({'Eigenvalues', 'Selected'}, 'Location', 'best', 'FontSize', 9)
title('Eigenvalue spectrum', 'Interpreter', 'latex', 'FontSize', 13)
grid on; box on;

nexttile
plot(abs(squeeze(StabRes.u(j,:,i))), StabGrid.y, 'k-', 'LineWidth', 1.5); hold on;
plot(abs(squeeze(StabRes.v(j,:,i))), StabGrid.y, 'b--', 'LineWidth', 1.2);
xlabel('$|\hat{q}|$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$\eta$', 'Interpreter', 'latex', 'FontSize', 12)
legend('$|\hat{u}|$', '$|\hat{v}|$', 'Interpreter', 'latex', 'Location', 'best')
title('Velocity eigenfunctions', 'Interpreter', 'latex', 'FontSize', 13)
grid on; box on;

nexttile
plot(abs(squeeze(StabRes.p(j,:,i))), StabGrid.y, 'g-', 'LineWidth', 1.5); hold on;
plot(abs(T_hat)/ampltd,               StabGrid.y, 'r-', 'LineWidth', 1.5);
plot(abs(rho_hat)/ampltd,             StabGrid.y, 'm-', 'LineWidth', 1.5);
legend({'$|\hat{p}|$', '$|\hat{T}|$', '$|\hat{\rho}|$'}, ...
       'Interpreter', 'latex', 'Location', 'best')
xlabel('$|\hat{q}|$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$\eta$', 'Interpreter', 'latex', 'FontSize', 12)
title('Thermodynamic eigenfunctions', 'Interpreter', 'latex', 'FontSize', 13)
grid on; box on;

end
