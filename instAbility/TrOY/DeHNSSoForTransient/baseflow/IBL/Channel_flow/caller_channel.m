clc; clear all; format long
restoredefaultpath 
%This part heere is just for saving%
Case_no=  102;
caseFolder = sprintf('Channel_Case_%03d', Case_no);
rootdir = fullfile(fileparts(mfilename('fullpath')), '..','..','..');
addpath(genpath(fullfile(rootdir, 'src')))





%% Inflow inputs %% 
% Note : In channel flow the velocity/length quantities are given already
% in a dimensional form% 

V = 1 ; 
nu  = 1E-03;  
nx = 1001 ;  L =120 ; 
H = 2;   S = 0 ; %(irrelevant the H,S in channel) 
N_cheb = 75; 
clear BF;
BF.model = 0;                             % choose between incompressible or compressible boundary layer (0: incompressible; 1: compressible)
% Mesh
BF.S     = S;                       %: (streamwise) domain start                            [m]       Scalar.
BF.L     = L;                        %: (streamwise) domain end                              [m]       Scalar.
BF.H     = H;                        %: (wall-normal) domain height                          [m]       Scalar.
BF.y_i   = 0;  %irrelevant for channel                        %: (wall-normal) Chebyshev node coordinate median       [m]       Scalar.                  
BF.nx    = nx;                        %: number of gridpoints in the streamwise direction     [-]       Scalar.
BF.ny    = N_cheb;                        %: number of gridpoints in the wall-normal direction    [-]       Scalar
BF.Ue=ones(1,BF.nx)*V;  % [m/s] streamwise velocity distribution (for ZPG, Ue = cst)
BF.We=0;                % [m/s] spanwise velocity (for unswept wings, We = 0)  

BF.nu = nu;   % [m^2/s] Kinematic viscocity


%% Plotting 
%Plotting features 2plt.mainFig   = 1;                     %choose whether not to (0) or to (1) plot after running a solver
%% Part B: boundary-layer fields (Simple Poisseuile) 
% Output: BL quantities. 


u = zeros(BF.ny,1);
v =  zeros(BF.ny,1); %Should be zero in Poiseuille flow 
w =  zeros(BF.ny,1); %Should be zero in Poiseuille flow
[yvec,DM] = chebdif(BF.ny,2); % differentiation matrices
y=yvec ; y_all = yvec; 
D1 = DM(:,:,1) ;
D2 = DM(:,:,2) ; 
u = 1-y.^2 ;

Du(:)=D1*u(:);
DDu(:)=D2*u(:);
Dv(:)=D1*v(:);
DDv(:)=D2*v(:);
Dw(:)=D1*w(:);
DDw(:)=D2*w(:);    

BL.u = u; BL.v = v; BL.w = w; BL.y = y;
BL.Du = Du; BL.DDu = DDu; BL.Dv = Dv; BL.DDv = DDv; BL.Dw = Dw; BL.DDw = DDw;


%load('C:\Users\LocalAdmin.TUD1001683\Desktop\Antonis\Load_files\BL_pois.mat','BL')
%nx =2000 ; 
ny = length(BL.y);
%BL_dot.X   = zeros(1,2000);
BF.X     = linspace(BF.S,BF.L,BF.nx);
BF.Y    = BL.y ;  %this if for malik 
BF.U    = repmat(BL.u,1,BF.nx) ;
BF.V    = repmat(BL.v,1,BF.nx) ; 
BF.W    = repmat(BL.w,1,BF.nx) ; 
BF.dyU  = repmat(BL.Du',1,BF.nx) ; 
BF.ddyU = repmat(BL.DDu',1,BF.nx) ; 
BF.dyV  = repmat(BL.Dv',1,BF.nx) ;
BF.ddyV = repmat(BL.DDv',1,BF.nx) ; 
BF.dyW  = repmat(BL.Dw',1,BF.nx);
BF.ddyW = repmat(BL.DDw',1,BF.nx) ;

BF.dxU  = zeros(BF.ny,BF.nx); 
BF.dxV  = zeros(BF.ny,BF.nx); 
BF.dxW  = zeros(BF.ny,BF.nx);
BF.Y_i =  1;





%% Global reference scales {ref.} & non-dimensionalization

%% Global reference scales {ref.} & non-dimensionalization

%------------ Reference scales ------------
BF.Uref=V;                       % [m/s] reference velocity
BF.lref=H/2; % [m] reference length scale (Blasius length-scale)
%BF.nu = BL.nu;                 % [m^2/s]
BF.Re = BF.Uref*BF.lref/BF.nu; % [-] Reference Reynolds nbr Re_0, based on the reference length scale (Re_0=sqrt{Re_x})

%------------ (Global) non-dimensionalization ------------
% Channel -> already normalised %
% BF.X=BF.X;   
% BF.Y=BF.Y;
% % velocities:
% BF.U=BF.U;
% BF.V=BF.V;
% BF.W=BF.W;
% % first x derivatives:
% BF.dxU  = BF.dudx;       
% BF.dxV  = BF.dvdx; 
% BF.dxW  = BF.dwdx; 
% % first y derivatives:
% BF.dyU  = BF.dyU;         
% BF.dyV  = BF.dyV; 
% BF.dyW  = BF.dyW; 
% % second y derivatives:
% BF.ddyU  = BF.ddyU;     
% BF.ddyV  = BF.ddyV; 
% BF.ddyW  = BF.ddyW; 
% 


%% Post processing %%
% No function for channel$



%% Save 

out_dir = fullfile('..','..','benchmark' , caseFolder , 'Clean','output_bf');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
%save(fullfile(out_dir, 'bf_clean.mat'), 'BF');
save(fullfile(out_dir,'bf_clean.mat'), 'BF');



% Checking % 

% contourf(BF.X,BF.Y,BF.ddyU)
% colorbar;