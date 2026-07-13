function [StabRes] = ILST(j,i,StabGrid,StabRes,bc_top,Opt)
%ILST  Compute inflow condition for mode j via local stability theory.
%
%   StabRes = ILST(j, i, StabGrid, StabRes)
%
%   Solves the Orr-Sommerfeld/Squire eigenvalue problem at station i,
%   filters the eigenvalue spectrum for the physical mode, normalises the
%   eigenfunction by max(|u|), and stores the result in StabRes.
%
%   Inputs:
%     j         — mode index
%     i         — streamwise station (typically 1 for inflow)
%     StabGrid  — grid structure (U, W, dyU, dyW, y, dzdy, d2zdy2, Re)
%     StabRes   — stability results (updated in-place with u, v, w, p, phi)

%% pre define some values from inputs

Re = StabGrid.Re;
U = StabGrid.U(:,i);
W = StabGrid.W(:,i);
dyU = StabGrid.dyU(:,i);
dyW = StabGrid.dyW(:,i);

omega = StabRes.omegavec(j);
beta = StabRes.betavec(j);
ny = StabGrid.ny;
dzdy = StabGrid.dzdy;
d2zdy2 = StabGrid.d2zdy2;

%% Find initial condition as solution to local stability eigenvalue problem

% Solve local eigenvalue problem ILST
if ~exist('bc_top','var'); bc_top = 'H_DR'; end
if ~exist('Opt','var');   Opt.plot = 'on'; end
if ~isfield(Opt,'verbose'); Opt.verbose = 'on'; end

alpha_shift = [];
if isfield(Opt,'alpha_shift') && ~isempty(Opt.alpha_shift)
    alpha_shift = Opt.alpha_shift;
end

if strcmpi(Opt.verbose,'on')
    fprintf(' Calculating IC for nonzero modes with ILST: ')
end
[eig,u,v,w,p] = ILST_solver(Re,flip(U),flip(W),flip(dyU),flip(dyW),omega,beta,ny,flip(dzdy),flip(d2zdy2),bc_top,alpha_shift);

% Find d99 for filtering purposes
k = find(U<=0.99*max(U),1,'first');
d99(i) = StabGrid.y(k);

% Select eigenvalue: EVfilter for full spectrum, direct for eigs shift
if isempty(alpha_shift)
    [alpha, index] = EVfilter(eig, v, StabGrid.y, omega, beta, d99);
else
    alpha = eig;   % scalar from eigs
    index = 1;
end

if strcmpi(Opt.verbose,'on')
    fprintf('alpha = %.6f %+.6fi\n', real(alpha), imag(alpha));
end

%% Normalize output by max(u)

StabRes.alpha(j,i) = alpha;

% Calculate new maximum amplitude after spline interpolation
% Select correct u
u                 = u(:,index);

% % Calculate amplitude
y_mdInt           = linspace(StabGrid.y(1),StabGrid.y(end),4000);
u1int             = interp1(StabGrid.y,abs(u),y_mdInt,'spline');
[ampltd, ~]       = max(abs(u1int));
 
% Normalize result by u_max
StabRes.u(j,:,i)  = u/ampltd;
StabRes.v(j,:,i)  = v(:,index)/ampltd;
StabRes.w(j,:,i)  = w(:,index)/ampltd;
StabRes.p(j,:,i)  = p(:,index)/ampltd;
 
% Create phi from eigenvectors
StabRes.phi(j,:,i) = ([ StabRes.u(j,:,i).';
                        StabRes.v(j,:,i).';
                        StabRes.w(j,:,i).';
                        StabRes.p(j,:,i).']);

%% Visual Check for Mode Selection

if ~strcmpi(Opt.plot, 'on'); return; end

figure(1); clf; set(gcf, 'Position', [200 200 1200 350]);
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% Base flow profiles at inflow
nexttile
plot(U, StabGrid.y, 'b-', 'LineWidth', 1.5); hold on;
plot(W, StabGrid.y, 'r-', 'LineWidth', 1.5);
plot(dyU, StabGrid.y, 'b--', 'LineWidth', 1.0);
plot(dyW, StabGrid.y, 'r--', 'LineWidth', 1.0);
xlabel('Velocity / derivative', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$\eta$', 'Interpreter', 'latex', 'FontSize', 12)
legend('U', 'W', 'dU/d\eta', 'dW/d\eta', 'Location', 'best', 'FontSize', 9)
title('Inflow base flow', 'Interpreter', 'latex', 'FontSize', 13)
grid on; box on;

% Eigenvalue spectrum
nexttile
plot(real(eig), imag(eig), 'k.', 'MarkerSize', 12); hold on;
plot(real(alpha), imag(alpha), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
alpha_r_center = omega;  % entropy branch at c=1, so alpha_r = omega
xlim([alpha_r_center - 0.1, alpha_r_center + 0.1])
ylim([-0.1, 0.02])
xlabel('$\alpha_r$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$\alpha_i$', 'Interpreter', 'latex', 'FontSize', 12)
legend({'Eigenvalues', 'Selected'}, 'Location', 'best', 'FontSize', 9)
title('Eigenvalue spectrum', 'Interpreter', 'latex', 'FontSize', 13)
grid on; box on;

% Shape function
nexttile
plot(abs(squeeze(StabRes.u(j,:,1))), StabGrid.y, 'k-', 'LineWidth', 1.5);
xlabel('$|\hat{u}|$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$\eta$', 'Interpreter', 'latex', 'FontSize', 12)
title('Eigenfunction', 'Interpreter', 'latex', 'FontSize', 13)
grid on; box on;

sgtitle(sprintf('Inflow ILST: $\\alpha = %.4f %+.4fi$', real(alpha), imag(alpha)), ...
    'Interpreter', 'latex', 'FontSize', 14);
drawnow;

end



