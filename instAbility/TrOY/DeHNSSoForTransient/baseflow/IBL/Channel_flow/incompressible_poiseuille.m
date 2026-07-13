function [BL,yvec,y_all,u, D1, D2] = incompressibleBL_poiseuille(BL)

%% Rename input variables as main structure

ny = BL.ny; H=BL.H ; y_i = BL.y_i; H = BL.H; nu = BL.nu; Ue = BL.Ue; 

%% Initialize fields

u = zeros(ny,1);
v = zeros(ny,1); %Should be zero in Poiseuille flow 
w = zeros(ny,1); %Should be zero in Poiseuille flow

% Initialize differentiation matrices
%[yvec,DM] = chebdif(ny,2); % differentiation matrices
%[y,D1,D2] = MappingMalik(H,y_i,yvec',DM(:,:,1),DM(:,:,2));

[yvec,DM] = chebdif(ny,2); % differentiation matrices
%[y,D1,D2] = MappingMalik(H,y_i,yvec',DM(:,:,1),DM(:,:,2));
%[y,D1,D2] = MappingMalik_pois(H,y_i,yvec',DM(:,:,1),DM(:,:,2));
%y = y/BL.H;

y=yvec ; y_all = yvec; 
D1 = DM(:,:,1) ;
D2 = DM(:,:,2) ; 
u = 1-y.^2 ;


% y_mirror = -y ;
% y_flip = flip(y_mirror);
% y_flip = y_flip(2:end);
% %reverse_y = -flip(y)
% y_all = [y; y_flip];
% uu = 1-y_all.^2;
% %y = reverse_y
% %clear DM yvec

%% Outputs

% Produce additional derivatives 
Du(:)=D1*u(:);
DDu(:)=D2*u(:);
Dv(:)=D1*v(:);
DDv(:)=D2*v(:);
Dw(:)=D1*w(:);
DDw(:)=D2*w(:);    


%% Rename output variables as main structure

BL.u = u; BL.v = v; BL.w = w; BL.y = y;
BL.Du = Du; BL.DDu = DDu; BL.Dv = Dv; BL.DDv = DDv; BL.Dw = Dw; BL.DDw = DDw;


end