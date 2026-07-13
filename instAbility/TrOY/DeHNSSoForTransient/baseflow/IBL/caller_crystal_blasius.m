
clc; clear all; format long
BL.model = 0;                             % choose between incompressible or compressible boundary layer (0: incompressible; 1: compressible)


%% add paths %% 



curr_dir = 'C:\Users\amrsynodinos\OneDrive - Delft University of Technology\Desktop\Codes';
curr_dir_ILST = fullfile(curr_dir,'repository-main-Numerical\Numerical');
curr_dir = 'C:\Users\amrsynodinos\OneDrive - Delft University of Technology\Desktop\Codes';
curr_dir_HNS = fullfile(curr_dir,'DeHNSSo NPSE Codes Sven Final');

%%Load some paths 07022025


% addpath((curr_dir_ILST));
% addpath(fullfile(curr_dir_ILST,'PartII_BaseFlow'));
% addpath(fullfile(curr_dir_ILST,'PartIII_Stability'));
% addpath(fullfile(curr_dir_ILST,'PartIII_Stability\Auxiliary'));
% addpath(fullfile(curr_dir_ILST,'PartII_BaseFlow\Similarity'));
% addpath(fullfile(curr_dir_ILST,'Mathematics'));
% addpath(fullfile(curr_dir_ILST,'PartII_BaseFlow\Auxiliary'));


addpath(fullfile(curr_dir_HNS,'Callers'));
addpath(fullfile(curr_dir_HNS,'Callers\TS-Waves in a Blasius BL'));
% addpath(fullfile(curr_dir_HNS,'Tools\Stability\HNS'));
% addpath(fullfile(curr_dir_HNS,'Tools\Stability\Nonlinear Prerequisites'));
% addpath(fullfile(curr_dir_HNS,'Tools\Grid generation'));
% addpath(fullfile(curr_dir_HNS,'Tools\Mathematics'));



%% y-parameteres, Height of the blob parameters %% 
H = 0.1;  
node_blob_c = 40; %like in ILST-OS %205 Ncheb=100;
N_cheb = 75;
h_target = 6; 

yi_factor = linspace(1,10000,10000); %50
V = 20.9 ; S = 0.1; nu = 1.51875e-5;
nx = 4001 ;  L =1.4 ; 
yi_target = Find_yi_h_target(h_target,node_blob_c,H,N_cheb,S,nu,V);
nodes_blob_all = node_blob_c; %28/04/2026

%% Crystal parameters 
lamda_y = 1;
Ampl = 0.01;
lcr_all = [0.01,0.02,0.05,0.1,0.15]; 
Smooth_waves_transition_all = 0; 
First_Clean_or = 0.25;
Last_Clean_or  = 0.25 + 0.15*L;
end_cr = 0.5; %for opposite ending%
% Modify slightly lcr such that we have an integer number of points per uc. 
dx=(L-S)/(nx-1);
for cr = 1:length(lcr_all)
   lcr = lcr_all(cr);
   ratio  = lcr/dx;
   int_ratio = ceil(ratio); % round up (no matters)
   lcr_new = int_ratio*dx;
   lcr_all_new(cr) = lcr_new; 
end 
lcr_all = lcr_all_new;
Nx_lcr = lcr_all/dx;



%Output : BF_mod{crystal_node_counter} ; 




%% Quasi-parallel inputs %%
xloc_int =  false   ;  Crystals_x_manual = linspace(0.2,0.9,8)     ;
Smooth_waves_transition_subgrid = 0;
N_mod = 7 ;
N_up =  1 ;
N_ds = 3  ; 
N_total = N_mod + N_up + N_ds ;
step_cr =1 ;


%% Generic inputs

% Mesh
BL.S     = S;                       %: (streamwise) domain start                            [m]       Scalar.
BL.L     = L;                        %: (streamwise) domain end                              [m]       Scalar.
BL.H     = H;                        %: (wall-normal) domain height                          [m]       Scalar.
BL.y_i   = yi_target;                       %: (wall-normal) Chebyshev node coordinate median       [m]       Scalar.                  
BL.nx    = nx;                        %: number of gridpoints in the streamwise direction     [-]       Scalar.
BL.ny    = N_cheb;                        %: number of gridpoints in the wall-normal direction    [-]       Scalar


BL.Ue=ones(1,BL.nx)*V;  % [m/s] streamwise velocity distribution (for ZPG, Ue = cst)
BL.We=0;                % [m/s] spanwise velocity (for unswept wings, We = 0)  

BL.nu = nu;   % [m^2/s] Kinematic viscocity

%% Plotting

% Plotting features
plt.mainFig   = 1;                     %choose whether not to (0) or to (1) plot after running a solver

%% Part B: boundary-layer fields

%[BL] = incompressibleBL_B1(BL);
%Antonis;
[BF] = incompressibleBL_B1(BL);

%% Global reference scales {ref.} & non-dimensionalization

%------------ Reference scales ------------
BF.Uref=V;                       % [m/s] reference velocity
BF.lref=sqrt(BL.S*BL.nu/BF.Uref); % [m] reference length scale (Blasius length-scale)
BF.nu = BL.nu;                 % [m^2/s]
BF.Re = BF.Uref*BF.lref/BF.nu; % [-] Reference Reynolds nbr Re_0, based on the reference length scale (Re_0=sqrt{Re_x})

%------------ (Global) non-dimensionalization ------------
% lengthscales:
BF.X=BL.x/BF.lref;   
BF.Y=BL.y/BF.lref;
BF.Y_i = BL.y_i/BF.lref;
BF.Y_i = BL.y_i/BF.lref;

% velocities:
BF.U=BL.u/BF.Uref;
BF.V=BL.v/BF.Uref;
BF.W=BL.w/BF.Uref;
% first x derivatives:
BF.dxU  = BL.dudx*(BF.lref/BF.Uref);       
BF.dxV  = BL.dvdx*(BF.lref/BF.Uref); 
BF.dxW  = BL.dwdx*(BF.lref/BF.Uref); 
% first y derivatives:
BF.dyU  = BL.Du*(BF.lref/BF.Uref);         
BF.dyV  = BL.Dv*(BF.lref/BF.Uref); 
BF.dyW  = BL.Dw*(BF.lref/BF.Uref); 
% second y derivatives:
BF.ddyU  = BL.DDu*(BF.lref)^2/BF.Uref;     
BF.ddyV  = BL.DDv*(BF.lref)^2/BF.Uref; 
BF.ddyW  = BL.DDw*(BF.lref)^2/BF.Uref; 

%% Postprocessing

[BF] = post_BL(BF);

%% Ref - Global length scales %%

Ref.U0 = BF.Uref ; Ref.L0 =BF.lref ; Ref.T0 = 293;  Ref.nu = BF.nu ;  Ref.Re = BF.Re;
lcr_all_dim = lcr_all/Ref.L0;

[BF_mod,First_Clean_dim_all,Last_Clean_dim_all]  =  Define_BF_mod(S,nx,Ref,BF,Ampl,lamda_y,lcr_all,lcr_all_dim,nodes_blob_all,Smooth_waves_transition_all,L,First_Clean_or,Last_Clean_or,end_cr);




% Plotting for checking 
BL_dot = BF_mod{1,1};
contourf(BL_dot.X,BL_dot.Y,BL_dot.U,'LineStyle','none')
ylim([0 20])
contourf(BL_dot.X,BL_dot.Y,BL_dot.U_per,'LineStyle','none')
colorbar
ylim([0 20])
% 


%% Save

out_dir = fullfile('..', 'output');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
save(fullfile(out_dir, 'bf_crystal.mat'), 'BF_mod');



















%% Save

out_dir = fullfile('..', 'output');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
save(fullfile(out_dir, 'bf_blasius.mat'), 'BF');
