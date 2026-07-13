
clc; clear all; format long

BL.model = 0;                             % choose between incompressible or compressible boundary layer (0: incompressible; 1: compressible)

%% Generic inputs

% Mesh
BL.S     = 0.1;                       %: (streamwise) domain start                            [m]       Scalar.
BL.L     = 1.7213;                        %: (streamwise) domain end                              [m]       Scalar.
BL.H     = 0.0585;                        %: (wall-normal) domain height                          [m]       Scalar.
BL.y_i   = BL.H/20;                       %: (wall-normal) Chebyshev node coordinate median       [m]       Scalar.                  
BL.nx    = 2000;                        %: number of gridpoints in the streamwise direction     [-]       Scalar.
BL.ny    = 75;     


% V = 20.9 ; S = 0.1; nu = 1.51875e-5; nx = 4001 ;  L =1.4 ; H = 0.1; 
% N_cheb = 75; 

% 
BL.S     = 0.1;                       %: (streamwise) domain start                            [m]       Scalar.
BL.L     = 1.4;                        %: (streamwise) domain end                              [m]       Scalar.
BL.H     = 0.1;                        %: (wall-normal) domain height                          [m]       Scalar.
BL.y_i   = BL.H/40;                       %: (wall-normal) Chebyshev node coordinate median       [m]       Scalar.                  
BL.nx    = 4001;                        %: number of gridpoints in the streamwise direction     [-]       Scalar.
BL.ny    = 75;     



%: number of gridpoints in the wall-normal direction    [-]       Scalar

% BL.nx  = 4000 ;
% BL.ny  = 50   ;

%V = 10;   % [m/s] freestream velocity at inlet
V = 20.9;   % [m/s] freestream velocity at inlet
BL.Ue=ones(1,BL.nx)*V;  % [m/s] streamwise velocity distribution (for ZPG, Ue = cst)
BL.We=0;                % [m/s] spanwise velocity (for unswept wings, We = 0)  

BL.nu = 1.511e-5;   % [m^2/s] Kinematic viscocity

%% Plotting

% Plotting features
plt.mainFig   = 1;                     %choose whether not to (0) or to (1) plot after running a solver

%% Part B: boundary-layer fields

[BL] = incompressibleBL_B1(BL);

%% Global reference scales {ref.} & non-dimensionalization

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

[BL] = post_BL(BL);

%% Save

out_dir = fullfile('..', 'output', 'benchmark');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
save(fullfile(out_dir, 'bf_blasius.mat'), 'BF');
