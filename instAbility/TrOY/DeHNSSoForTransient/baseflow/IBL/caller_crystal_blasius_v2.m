
clc; clear all; format long
restoredefaultpath;

%% add paths %% 

Case_no=  1;
caseFolder = sprintf('Case_%03d', Case_no);




curr_dir = 'C:\Users\amrsynodinos\OneDrive - Delft University of Technology\Desktop\Codes';
curr_dir_ILST = fullfile(curr_dir,'repository-main-Numerical\Numerical');
curr_dir = 'C:\Users\amrsynodinos\OneDrive - Delft University of Technology\Desktop\Codes';
curr_dir_HNS = fullfile(curr_dir,'DeHNSSo NPSE Codes Sven Final');
addpath(fullfile(curr_dir_HNS,'Callers'));
addpath(fullfile(curr_dir_HNS,'Callers\TS-Waves in a Blasius BL'));


%% Inflow inputs %%

V = 20.9 ; S = 0.1; nu = 1.51875e-5; nx = 4001 ;  L =1.4 ; H = 0.1; 
N_cheb = 75; 

%% y-parameteres, Height of the blob parameters %% 

%% y-Crystal parameters %%
N_blob = 40;
lamda_y = 1;
Ampl =0.01;
Ampl_all = [0.01,0.02];
Ampl_all =[0.01];
yi_factor = linspace(1,10000,10000); %50
%node_blob_c = 40; %like in ILST-OS %205 Ncheb=100;
N_cheb = 75;
h_target_all = [0.5,1,6];
h_target_all = [1];

for j = 1:length(h_target_all)
    h_target = h_target_all(j);
    yi_target = Find_yi_h_target(h_target,N_blob,H,N_cheb,S,nu,V);
    nodes_blob_all(j) = N_blob   ; %28/04/2026%
    yi_target_all(j) = yi_target ;
end 

%% x-Crystal parameters 
lcr_all = [0.01,0.02,0.05,0.1,0.15]; 
lcr_all = [0.01,0.1,0.15];
lcr_all = [0.01];
Smooth_waves_transition_all = 0; 
First_Clean_or = 0.25;
Last_Clean_or  = 0.25 + 0.15*L;
end_cr = 0.5; %for opposite ending
end_cr = 0;
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


%% Iterations for BF, because the yi is modified 

for j = 1:length(yi_target_all)
    clear BF;
    BF.model = 0;                             % choose between incompressible or compressible boundary layer (0: incompressible; 1: compressible)
    yi_target = yi_target_all(j);
% Mesh
    BF.S     = S;                       %: (streamwise) domain start                            [m]       Scalar.
    BF.L     = L;                        %: (streamwise) domain end                              [m]       Scalar.
    BF.H     = H;                        %: (wall-normal) domain height                          [m]       Scalar.
    BF.y_i   = yi_target;                       %: (wall-normal) Chebyshev node coordinate median       [m]       Scalar.                  
    BF.nx    = nx;                        %: number of gridpoints in the streamwise direction     [-]       Scalar.
    BF.ny    = N_cheb;                        %: number of gridpoints in the wall-normal direction    [-]       Scalar


    BF.Ue=ones(1,BF.nx)*V;  % [m/s] streamwise velocity distribution (for ZPG, Ue = cst)
    BF.We=0;                % [m/s] spanwise velocity (for unswept wings, We = 0)  

    BF.nu = nu;   % [m^2/s] Kinematic viscocity
    
    %% Plotting
    
    % Plotting features
    plt.mainFig   = 1;                     %choose whether not to (0) or to (1) plot after running a solver
    
    %% Part B: boundary-layer fields
    
    %[BF] = incompressibleBL_B1(BF);
    %Antonis;
    [BF] = incompressibleBL_B1(BF);
    
    %% Global reference scales {ref.} & non-dimensionalization
    
    %------------ Reference scales ------------
    BF.Uref=V;                       % [m/s] reference velocity
    BF.lref=sqrt(BF.S*BF.nu/BF.Uref); % [m] reference length scale (Blasius length-scale)
    BF.nu = BF.nu;                 % [m^2/s]
    BF.Re = BF.Uref*BF.lref/BF.nu; % [-] Reference Reynolds nbr Re_0, based on the reference length scale (Re_0=sqrt{Re_x})
    
    %------------ (Global) non-dimensionalization ------------
    % lengthscales:
    BF.X=BF.x/BF.lref;   
    BF.Y=BF.y/BF.lref;
    BF.Y_i = BF.y_i/BF.lref;
    BF.Y_i = BF.y_i/BF.lref;
    
    % velocities:
    BF.U=BF.u/BF.Uref;
    BF.V=BF.v/BF.Uref;
    BF.W=BF.w/BF.Uref;
    % first x derivatives:
    BF.dxU  = BF.dudx*(BF.lref/BF.Uref);       
    BF.dxV  = BF.dvdx*(BF.lref/BF.Uref); 
    BF.dxW  = BF.dwdx*(BF.lref/BF.Uref); 
    % first y derivatives:
    BF.dyU  = BF.Du*(BF.lref/BF.Uref);         
    BF.dyV  = BF.Dv*(BF.lref/BF.Uref); 
    BF.dyW  = BF.Dw*(BF.lref/BF.Uref); 
    % second y derivatives:
    BF.ddyU  = BF.DDu*(BF.lref)^2/BF.Uref;     
    BF.ddyV  = BF.DDv*(BF.lref)^2/BF.Uref; 
    BF.ddyW  = BF.DDw*(BF.lref)^2/BF.Uref; 
    
    %% Postprocessing
    
    [BF] = post_BL(BF);
    
    %% Ref - Global length scales %%
    
    Ref.U0 = BF.Uref ; Ref.L0 =BF.lref ; Ref.T0 = 293;  Ref.nu = BF.nu ;  Ref.Re = BF.Re; %(Could be outside iteration, but works and here) 
    lcr_all_dim = lcr_all/Ref.L0;
    for k = 1:length(Ampl_all);
        Ampl = Ampl_all(k);
        for i = 1:length(lcr_all)
            %[BF_mod,First_Clean_dim_all,Last_Clean_dim_all]  =  Define_BF_mod(S,nx,Ref,BF,Ampl,lamda_y,lcr_all,lcr_all_dim,nodes_blob_all,Smooth_waves_transition_all,L,First_Clean_or,Last_Clean_or,end_cr);
            [BF_mod,First_Clean_dim,Last_Clean_dim]  =  Define_BF_mod_v2(S,nx,Ref,BF,Ampl,lamda_y,lcr_all,nodes_blob_all,Smooth_waves_transition_all,L,First_Clean_or,Last_Clean_or,end_cr,i,j);
            
            BF_mod_all{k,j,i} = BF_mod;
            First_clean_dim_all(i) = First_Clean_dim;
            Last_clean_dim_all(i) = Last_Clean_dim;
            
        end 
    end %k;   
end 


%The organisation is : Ampl-blob-wavelength %

%%
% Plotting for checking 
% BL_dot = BF_mod_all{2,2,3};
% contourf(BL_dot.X,BL_dot.Y,BL_dot.U,'LineStyle','none')
% ylim([0 20])
% contourf(BL_dot.X,BL_dot.Y,BL_dot.U_per,'LineStyle','none')
% colorbar
% ylim([0 20])
% % 


%% Save

%out_dir = fullfile('..', 'output');
out_dir = fullfile('..','benchmark/',caseFolder);

if ~exist(out_dir, 'dir'); mkdir(out_dir); end
% %save(fullfile(out_dir,'benchmark/' ,'bf_crystal.mat'), 'BF_mod_all'); % maybe not necessary%
% save(fullfile(out_dir,'benchmark/' ,'Ampl_all.mat'), 'Ampl_all');
% save(fullfile(out_dir,'benchmark/' ,'h_target_all.mat'), 'h_target_all');
% save(fullfile(out_dir,'benchmark/' ,'lcr_all.mat'), 'lcr_all');

save(fullfile(out_dir,'Ampl_all.mat'), 'Ampl_all');
save(fullfile(out_dir,'h_target_all.mat'), 'h_target_all');
save(fullfile(out_dir,'lcr_all.mat'), 'lcr_all');

% One more step, create folders with the specific cases and save them
% inside there.

for k = 1:length (Ampl_all)
  for j = 1:length(h_target_all) 
      for i = 1:length(lcr_all)
         Ampl = Ampl_all(k); h_target = h_target_all(j); lcr = lcr_all(i);
         foldername =sprintf('A%.2f_h%.2f_lcr%.2f', Ampl ,h_target,lcr);
         %folder_path = fullfile(out_dir,'benchmark/',foldername);
        folder_path = fullfile(out_dir,foldername,'output_bf');

         if ~exist(folder_path, 'dir')
            mkdir(folder_path)
         end
         bf_crystal = BF_mod_all{k,j,i};
         save( fullfile(folder_path, 'bf_crystal.mat'), 'BF');
         save( fullfile(folder_path, 'Ampl.mat'), 'Ampl');
         save( fullfile(folder_path, 'h_target.mat'), 'h_target');
         save( fullfile(folder_path, 'lcr.mat'), 'lcr');

         


      end
  end   
end 













