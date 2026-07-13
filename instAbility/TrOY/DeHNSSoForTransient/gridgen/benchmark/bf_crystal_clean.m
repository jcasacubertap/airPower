%% Gridgen benchmark — Blasius BL base flow
%  Cartesian structured input, already at DeHNSSo grid resolution.
%  Passed through unchanged (no resampling).

clear all; close all; clc
restoredefaultpath
rootdir = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(rootdir, genpath(fullfile(rootdir, 'src')))

%% Input

Case_no=  2;
caseFolder = sprintf('Case_%03d', Case_no);



%% 

input.folder_case     = fullfile(rootdir, '..', 'baseflow', 'benchmark',caseFolder);
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
            input.folder     = fullfile(input.folder_case,foldername,'output_bf');
            input.filename   = 'bf_crystal.mat';
            input.format     = 'mat';
            input.structured = true;


                    %% Resampling parameters
            % Pass-through: leave n_eta_new empty so the Cartesian branch keeps the
            % source grid as-is (bf_blasius already matches DeHNSSo's target).
            params.n_eta_new = [];
            params.n_xi_new  = [];
            params.rescale   = false;    % already non-dim
            params.plot      = false;
            
            %% Grid-generation options
            params.FD1_order    = 4;
            params.FD2_order    = 2;
            params.eta_method   = 'cheb';
            params.FD_eta_order = 4;
            %params.plot =1 ; Antonis for plotting. 
            
            %% Output
            %output.folder   = fullfile(rootdir, '..', 'input');
            output.folder   = fullfile( input.folder, '..', 'input_stab');
            %output.folder   = fullfile( input.folder, '..', 'input_stab');
            if ~exist( output.folder, 'dir') 
                mkdir(folder_path) 
            end
            output.filename = 'DeHNSSo_input.mat';
            fprintf('DeHNSSo_input created: %s\n', foldername)

            
            %% Run
            StabGrid = main_gridgen(input, params, output);



        end %i
    end %j
end %k


% Plotting for checking 
% BL_dot = BF;
% contourf(StabGrid.x,StabGrid.y,StabGrid.dyU,'LineStyle','none')
% ylim([0 20])
% colorbar;
% contourf(BL_dot.X,BL_dot.Y,BL_dot.V,'LineStyle','none')
% colorbar
% ylim([0 20])





% Do i need to consider more inputs here? Like for example No_Ampl, No_lcr,
% and No_blobs




