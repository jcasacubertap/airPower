
clc; clear all; format long
restoredefaultpath
Case_no=  2;
caseFolder = sprintf('Case_%03d', Case_no);



%% Inflow inputs %%

V = 20.9 ; S = 0.1; nu = 1.511e-5; nx = 4001 ;  L =1.4 ; H = 0.1; 
N_cheb = 75; 
clear BF;
BF.model = 0;                             % choose between incompressible or compressible boundary layer (0: incompressible; 1: compressible)
% Mesh
BF.S     = S;                       %: (streamwise) domain start                            [m]       Scalar.
BF.L     = L;                        %: (streamwise) domain end                              [m]       Scalar.
BF.H     = H;                        %: (wall-normal) domain height                          [m]       Scalar.
BF.y_i   = BF.H/40;                       %: (wall-normal) Chebyshev node coordinate median       [m]       Scalar.                  
BF.nx    = nx;                        %: number of gridpoints in the streamwise direction     [-]       Scalar.
BF.ny    = N_cheb;                        %: number of gridpoints in the wall-normal direction    [-]       Scalar



% 
% BF.S     = 0.1;                       %: (streamwise) domain start                            [m]       Scalar.
% BF.L     = 1.7213;                        %: (streamwise) domain end                              [m]       Scalar.
% BF.H     = 0.0585;                        %: (wall-normal) domain height                          [m]       Scalar.
% BF.y_i   = BF.H/20;                       %: (wall-normal) Chebyshev node coordinate median       [m]       Scalar.                  
% BF.nx    = 2000;                        %: number of gridpoints in the streamwise direction     [-]       Scalar.
% BF.ny    = 75;     




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



%%
% Plotting for checking 
BL_dot = BF;
contourf(BL_dot.X,BL_dot.Y,BL_dot.dyU,'LineStyle','none')
ylim([0 20])
contourf(BL_dot.X,BL_dot.Y,BL_dot.V,'LineStyle','none')
colorbar
ylim([0 20])
% 

%%
%% Save

out_dir = fullfile('..','benchmark' , caseFolder , 'Clean','output_bf');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
%save(fullfile(out_dir, 'bf_clean.mat'), 'BF');
save(fullfile(out_dir,'bf_clean.mat'), 'BF');





















%% Save

out_dir = fullfile('..', 'output');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
save(fullfile(out_dir, 'bf_blasius.mat'), 'BF');
