%% DeHNSSo — Blasius boundary layer (TS waves)
%  Base flow: bf_blasius.mat (see gridgen/benchmark/blasius.m)

clear all; close all; clc
restoredefaultpath
rootdir = fullfile(fileparts(mfilename('fullpath')), '..','..');
%% Stability parameters 
frequency = 0;
beta      = 1;
%
Case_no=  102;
caseFolder = sprintf('Channel_Case_%03d', Case_no);


%% 

addpath(genpath(fullfile(rootdir, 'src')))
input.folder_case     = fullfile(rootdir,  'baseflow', 'benchmark',caseFolder,'Clean');
input.folder     = fullfile(input.folder_case,'input_stab');
load(fullfile(input.folder, 'DeHNSSo_input.mat'));   % provides StabGrid
StabGrid.case = 'Channel';
% StabGrid.dzdy = zeros(StabGrid.ny,1);
% StabGrid.d2zdy2  = zeros(StabGrid.ny,1);
% =========================================================================
%% Spectral content
% =========================================================================
%  total modes: (2M+1)(2N+1), including conjugates

Stab.M = 1;                 % number of omega (temporal) modes
Stab.N = 0;                 % number of beta (spanwise) modes
% Stab.omega_0 = 0.034576;    % fundamental angular frequency
%Stab.omega_0 = 2*pi*frequency*(StabGrid.lref/StabGrid.Uref);    % fundamental angular frequency
Stab.omega_0 = frequency ;    % already dimensionalised. 
Stab.beta_0  = beta;           % fundamental spanwise wavenumber (0 = 2D)

% =========================================================================
%% Inflow conditions
% =========================================================================

%Stab.IC       = 'ILST';     % 'ILST' (local stability) or 'LOAD' (from file)
Stab.phaseRef = 'pwall';    % 'pwall' (TS) or 'umax' (CFI)

Stab.IC       = 'LOAD'; 
Stab.ICfile = fullfile(input.folder, 'StabRes_previous.mat');  % required for LOAD

Opt.linear   = 'on';        % 'on' linear (NLT off, 1 iter) | 'off' nonlinear
Stab.A0_fund = 0.00125*sqrt(2);   % fundamental disturbance amplitude

% =========================================================================
%% Boundary conditions
% =========================================================================
%
% Wall (y=0): always inhomogeneous Dirichlet.
%   Stab.bcw(k,ix,j) — wall BC value.  k=1 u, k=2 v, k=3 w.
%   Default: zero (no-slip).
%
% Freestream (y=H): per-component BC type.
%   Opt.bc_top{k} — type for component k=1 u, k=2 v, k=3 w.
%     'H_DR'  = homogeneous Dirichlet   (q = 0, default)
%     'IH_DR' = inhomogeneous Dirichlet (q = bct value, requires Stab.bct)
%     'H_NM'  = homogeneous Neumann     (dq/dy = 0)
%   Pressure: no BC at freestream (determined by the equations).

Opt.bc_top = {'H_DR', 'H_DR', 'H_DR'};

% =========================================================================
%% Buffer — outflow treatment
% =========================================================================
%   'on'   = amplitude damping (scales solution -> 0 near outflow)
%   'para' = parabolisation (removes d2/dx2 in buffer -> no upstream reflections)
%   'off'  = no buffer

Opt.buffer = 'para';
Opt.xb     = 85;            % buffer start [% of domain]

% =========================================================================
%% Solver options
% =========================================================================

Opt.plot   = 'on';          % 'on' | 'off'
% Opt.Conv = 1e-4;          % final convergence criterion
% Opt.TH   = 1e-11;         % NLT activation threshold

% =========================================================================
%% Performance
% =========================================================================

Opt.lu_mode = 'full';       % LU cache: 'full' | 'auto' (LRU) | 'none'
Opt.parfor  = 'off';        % parallel assembly/solve: 'on' | 'off'
Opt.gpu     = 'off';        % GPU sparse solve: 'on' | 'off'

% =========================================================================
%% Nonlinear solver choice
% =========================================================================
% 'picard' — Picard with AF ramping + mode activation (default)
% 'newton' — JFNK with per-mode preconditioner

Opt.solver = 'picard';      % 'picard' (+ Newton fallback) | 'newton'
% Opt.newton_tol = 1e-10;

% =========================================================================
%% Run
% =========================================================================
% Opt.rerun = 'on' reuses LHS from a previous run (skips matrix assembly).

Opt.rerun = 'off';

tic_total = tic;
if strcmpi(Opt.rerun, 'on')
    [StabRes, StabGrid, LHS] = main(Stab, StabGrid, Opt, LHS);
else
    [StabRes, StabGrid, LHS] = main(Stab, StabGrid, Opt);
end
t_total = toc(tic_total);
w = whos;
fprintf('\n========================================\n');
fprintf('  Total execution time: %.2f s\n', t_total);
fprintf('  Workspace memory:     %.1f MB\n', sum([w.bytes])/1e6);
fprintf('========================================\n');

% =========================================================================
%% Post-processing
% =========================================================================
close all





output.folder = fullfile(input.folder,'..','output_stab');

if ~exist(output.folder, 'dir')
mkdir(output.folder)
end

output.Gridname_mat  = 'StabGrid.mat';
output.Resname_mat  =  'StabRes.mat';
output.Gridname = 'StabGrid';
output.Resname  =  'StabRes';
save(fullfile( output.folder, output.Gridname_mat), output.Gridname);
save(fullfile( output.folder, output.Resname_mat),  output.Resname );














u =squeeze(StabRes.u(2,:,:));
v =squeeze(StabRes.v(2,:,:));
w =squeeze(StabRes.w(2,:,:));
p =squeeze(StabRes.p(2,:,:));



x = StabGrid.x;
y = StabGrid.y;
A_u = max(abs(u));
A_v = max(abs(v));
A_p = max(abs(p));

%----------------- Real ------------------%
fig = figure;
fig.Position = [100 100 800 600];
subplot(4,1,1)
contourf(x,y,real(u)./A_u,'LineStyle','none')
colorbar;
title('$Real(u)$','Interpreter','latex')
set(gca,'FontSize',12)
subplot(4,1,2) 
contourf(x,y,real(v)./A_v,'LineStyle','none')
colorbar;
title('$Real(v)$','Interpreter','latex')

set(gca,'FontSize',12)
subplot(4,1,3) 
contourf(x,y,real(p)./A_p,'LineStyle','none')
title('$Real(p)$','Interpreter','latex')

colorbar;
set(gca,'FontSize',12)
subplot(4,1,4) 
plot(x,log(A_u/A_u(2)),'--k','Displayname','$N_{u}$')
hold on      
plot(x,log(A_v/A_v(1)),'--r','Displayname','$N_{v}$')
hold on      
%plot(x,log(A_p/A_p(1)),'--b','Displayname','$N_{p}$')
legend('Interpreter','latex')
ylabel('$N$','Interpreter','latex')
set(gca,'FontSize',12)

%----------------- Abs ------------------%

fig = figure;
fig.Position = [100 100 800 600];
subplot(4,1,1)
contourf(x,y,abs(u)./A_u,'LineStyle','none')
colorbar;
title('$Abs(u)$','Interpreter','latex')
set(gca,'FontSize',12)
subplot(4,1,2) 
contourf(x,y,abs(v)./A_v,'LineStyle','none')
colorbar;
title('$Abs(v)$','Interpreter','latex')

set(gca,'FontSize',12)
subplot(4,1,3) 
contourf(x,y,abs(p)./A_p,'LineStyle','none')
title('$Abs(p)$','Interpreter','latex')

colorbar;
set(gca,'FontSize',12)
subplot(4,1,4) 

plot(x,log(A_u/A_u(2)),'--k','Displayname','$N_{u}$')
hold on      
plot(x,log(A_v/A_v(1)),'--r','Displayname','$N_{v}$')
hold on      
%plot(x,A_p,'--b','Displayname','$N_{p}$')
legend('Interpreter','latex')
ylabel('$N$','Interpreter','latex')
set(gca,'FontSize',12)
% 
% plot(x,log(A_u/A_u(2)),'--k','Displayname','$N_{u}$')
% hold on      
% plot(x,log(A_v/A_v(1)),'--r','Displayname','$N_{v}$')
% hold on      
% plot(x,log(A_p/A_p(1)),'--b','Displayname','$N_{p}$')
% legend('Interpreter','latex')
% ylabel('$N$','Interpreter','latex')
% set(gca,'FontSize',12)


StabRes_v2 = StabRes;
StabGrid_v2 = StabGrid;
% BF_v2 = BF;
% Grid_v1 = Grid;
% Stab_v1 = Stab;
% Opt_v1 = Opt;

% Comparison_folder = 'C:\Users\amrsynodinos\OneDrive - Delft University of Technology\Desktop\DeHNSSo_new_valid';
% folder_v2 = fullfile(Comparison_folder,'v2');
% save(fullfile(folder_v2,'StabRes_v2.mat'),'StabRes_v2');
% save(fullfile(folder_v2,'StabGrid_v2.mat'),'StabGrid_v2');
% save(fullfile(folder_v1,'BF_v1.mat'),'BF_v1');
% save(fullfile(folder_v1,'Stab_v1.mat'),'Stab_v1');
% save(fullfile(folder_v1,'Opt_v1.mat'),'Opt_v1');

% Reference for overlay (set ref_file = [] to skip)
ref_file  = fullfile(rootdir, 'literature', 'StabRes_Blasius_NPSE.mat');
ref_label = 'NPSE';

if strcmpi(Opt.plot, 'on')
    if ~isempty(ref_file) && isfile(ref_file)
        plot_amplitudes(StabGrid, StabRes, struct('file', ref_file, 'label', ref_label));
    else
        if ~isempty(ref_file)
            fprintf('  [ref overlay skipped — file not found: %s]\n', ref_file);
        end
        plot_amplitudes(StabGrid, StabRes);
    end
    plot_stability(StabGrid, StabRes, Opt);
end

