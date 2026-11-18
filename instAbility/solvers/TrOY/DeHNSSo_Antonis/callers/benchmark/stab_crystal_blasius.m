%% DeHNSSo — Blasius boundary layer (TS waves)
%  Base flow: bf_blasius.mat (see gridgen/benchmark/blasius.m)

clear all; close all; clc
restoredefaultpath
rootdir = fullfile(fileparts(mfilename('fullpath')), '..','..');
%% Stability parameters 
frequency =285;
%
Case_no=  1;
caseFolder = sprintf('Case_%03d', Case_no);

%% 

addpath(genpath(fullfile(rootdir, 'src')))
input.folder_case     = fullfile(rootdir,  'baseflow', 'benchmark',caseFolder);
files = dir(input.folder_case );
dirFlags = [files.isdir];
folders = files(dirFlags);
folders = folders(~ismember({folders.name}, {'.', '..'}));
load (fullfile(input.folder_case ,'Ampl_all.mat'),'Ampl_all')
load (fullfile(input.folder_case,'h_target_all.mat'),'h_target_all')
load (fullfile(input.folder_case,'lcr_all.mat'),'lcr_all')

for k = 1:length(Ampl_all)
    for j = 1:length(h_target_all)
        for i = 1:length(lcr_all)
            Ampl    = Ampl_all(k);
            h_target = h_target_all(j);
            lcr      = lcr_all(i);
            foldername = sprintf('A%.2f_h%.2f_lcr%.2f', ...
                                  Ampl, h_target, lcr);
            folder_path = fullfile(input.folder_case, foldername);
            if exist(folder_path, 'dir')
                fprintf('EXISTS: %s\n', foldername)
            else
                fprintf('MISSING: %s\n', foldername)
            end
            input.folder     = fullfile(input.folder_case,foldername,'input_stab');
            load(fullfile(input.folder, 'DeHNSSo_input.mat'));   % provides StabGrid
        
            % =========================================================================
            %% Spectral content
            % =========================================================================
            %  total modes: (2M+1)(2N+1), including conjugates
            
            Stab.M = 1;                 % number of omega (temporal) modes
            Stab.N = 0;                 % number of beta (spanwise) modes
           % Stab.omega_0 = 0.034576;    % fundamental angular frequency
            Stab.omega_0 = 2*pi*frequency*(StabGrid.lref/StabGrid.Uref);    % fundamental angular frequency
            Stab.beta_0  = 0;           % fundamental spanwise wavenumber (0 = 2D)
            
            % =========================================================================
            %% Inflow conditions
            % =========================================================================
            
            Stab.IC       = 'ILST';     % 'ILST' (local stability) or 'LOAD' (from file)
            Stab.phaseRef = 'pwall';    % 'pwall' (TS) or 'umax' (CFI)
            % Stab.ICfile = fullfile(input.folder, 'StabRes_previous.mat');  % required for LOAD
            
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
            
            
            u =squeeze(StabRes.u(2,:,:));
            A_u = max(abs(u));
                  
            StabRes_v2 = StabRes;
            StabGrid_v2 = StabGrid;
            % BF_v2 = BF;
            % Grid_v1 = Grid;
            % Stab_v1 = Stab;
            % Opt_v1 = Opt;
            
            Comparison_folder = 'C:\Users\amrsynodinos\OneDrive - Delft University of Technology\Desktop\DeHNSSo_new_valid';
            folder_v2 = fullfile(Comparison_folder,'v2');
            save(fullfile(folder_v2,'StabRes_v2.mat'),'StabRes_v2');
            save(fullfile(folder_v2,'StabGrid_v2.mat'),'StabGrid_v2');
            % save(fullfile(folder_v1,'BF_v1.mat'),'BF_v1');
            % save(fullfile(folder_v1,'Stab_v1.mat'),'Stab_v1');
            % save(fullfile(folder_v1,'Opt_v1.mat'),'Opt_v1');
            end 
    end 
end
