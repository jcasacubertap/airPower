
clc; clear all; format long
caller_dir = fileparts(mfilename('fullpath'));
rootdir    = fullfile(caller_dir, '..', '..');     % refactor repo root
addpath(genpath(rootdir));

BL.model = 0;                             % choose between incompressible or compressible boundary layer (0: incompressible; 1: compressible)

%% Generic inputs

% Mesh
BL.S     = 0.2418;                       %: (streamwise) domain start                            [m]       Scalar.
BL.L     = 1.7213;                        %: (streamwise) domain end                              [m]       Scalar.
BL.H     = 0.0585;                        %: (wall-normal) domain height                          [m]       Scalar.
BL.y_i   = BL.H/20;                       %: (wall-normal) Chebyshev node coordinate median       [m]       Scalar.                  
BL.nx    = 800;                        %: number of gridpoints in the streamwise direction     [-]       Scalar.
BL.ny    = 40;                        %: number of gridpoints in the wall-normal direction    [-]       Scalar

V = 10;   % [m/s] freestream velocity at inlet
BL.Ue=ones(1,BL.nx)*V;  % [m/s] streamwise velocity distribution (for ZPG, Ue = cst)
BL.We=0;                % [m/s] spanwise velocity (for unswept wings, We = 0)  

BL.nu = 1.511e-5;   % [m^2/s] Kinematic viscocity

%% Plotting

% Plotting features
plt.mainFig   = 1;                     %choose whether not to (0) or to (1) plot after running a solver

%% Part B: boundary-layer fields

[BL] = incompressibleBL_B1(BL);

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
% (Derivatives are produced by gridgen — see Interpol/src/compute_bf_derivatives.m.)

%% Postprocessing

[BL] = post_BL(BL);

%% Save

outdir = fullfile(caller_dir, '..', 'output');
if ~exist(outdir, 'dir'); mkdir(outdir); end
save(fullfile(outdir, 'BF_Blasius.mat'), 'BF');
