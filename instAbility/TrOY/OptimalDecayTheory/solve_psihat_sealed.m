%% solve_psihat_mode.m
% Solves the single spanwise-Fourier-mode
% boundary-value problem for the cross-plane streamfunction,
%
%        psihat'' - beta^2 psihat = i beta (duB/dy) uhat,   psihat(+-1) = 0,
%
% for a prescribed streak mode  u'(y,z) = Re[ uhat(y) e^{i beta z} ],
% uhat(y) = A(y) e^{i phi(y)}].  The cross-stream velocity then
% follows as  vhat = i beta psihat,  what = -d psihat/dy.
%
% Discretisation: Chebyshev collocation in y on [-1,1]. A single
% wavenumber beta is solved. A general spanwise-periodic streak is recovered
% by calling this for each constituent beta and superposing.
%
% Base flow (default): plane Poiseuille  uB = 1 - y^2,  duB/dy = -2y.

clear; clc;

%% ----- inputs ------------------------------------------------------------
beta = 2.0;                 % spanwise wavenumber  
N    = 100;                 % Chebyshev resolution in y

%Libary Afun: @(y) sin(pi*y); @(y) (1-y.^2).*exp(-((y+0.7)/0.25).^2);
%@(y) (1-y.^2).*cos(2*pi*y); @(y) cos(pi*y/2); @(y) 1 - y.^2;           
Afun   = @(y) 1 - y.^2;     % streak amplitude  A(y)
phifun = @(y) 0.8*y;          % streak phase      phi(y)   (e.g. 0.8*y to tilt)
dUBdy  = @(y) -2*y;         % base-flow shear   duB/dy

%% ----- discretisation and assembly  ----------------------------
[D,y] = cheb(N);            % nodes y (col, +1..-1) and 1st-derivative matrix
D2    = D*D;
Up    = dUBdy(y);
uhat  = Afun(y).*exp(1i*phifun(y));          % uhat(y) = A e^{i phi}  (eq. 1.10)

L = D2 - beta^2*eye(N+1);                    % operator  d^2/dy^2 - beta^2
r = 1i*beta*(Up.*uhat);                      % rhs  i beta (duB/dy) uhat

% Dirichlet boundary conditions  psihat(+-1) = 0
L(1,:)   = 0; L(1,1)     = 1; r(1)   = 0;    % y = +1
L(N+1,:) = 0; L(N+1,N+1) = 1; r(N+1) = 0;    % y = -1

psihat = L\r;                                % solve main equation

%% ----- cross-stream velocity ----------------------------------
vhat =  1i*beta*psihat;                      % vhat = i beta psihat
what = -D*psihat;                            % what = -d psihat/dy

%% ----- verification ------------------------------------------------------
cont   = D*vhat + 1i*beta*what;              % continuity  (machine zero)
omx    = D*what - 1i*beta*vhat;              % streamwise vorticity, modal
omx_id = omx + Up.*(1i*beta*uhat);           % omega_x = -(duB/dy) du'/dz  -> ~0
fprintf('max |continuity|            : %.2e\n', max(abs(cont)));
fprintf('max |omega_x + U'' du''/dz|   : %.2e\n', max(abs(omx_id(2:N))));
fprintf('psihat at walls |+1|,|-1|   : %.2e , %.2e\n', abs(psihat(1)), abs(psihat(N+1)));

%% ----- reconstruct physical fields over two wavelengths -----------------
Lz = 2*pi/beta; Nz = 240; z = linspace(0,2*Lz,Nz);
E  = exp(1i*beta*z);
rec = @(a) real(a*E);
up  = rec(uhat);  psi = rec(psihat);  vp = rec(vhat);  wp = rec(what);

%% ----- plot: streak (colour) + roller streamlines + (w',v') quiver ------
[Z,Y] = meshgrid(z,y);
figure('Color','w','Position',[100 100 820 420]);
pcolor(Z,Y,up); shading interp; hold on
colormap(redblue_local()); cmax=max(abs(up(:))); caxis([-cmax cmax]);
cb=colorbar; cb.Label.String='u'' (streak)';
contour(Z,Y,psi,12,'k','LineWidth',0.5);
iy=1:max(1,round(N/16)):N+1; iz=1:max(1,round(Nz/18)):Nz;
quiver(Z(iy,iz),Y(iy,iz),wp(iy,iz),vp(iy,iz),0.9,'k');
xlabel('z'); ylabel('y'); axis tight
title(sprintf('Mode \\beta = %.2f : streak and de-streaking roller',beta));
hold off

%% ======================= local functions =================================
function [D,x] = cheb(N)
    if N==0, D=0; x=1; return; end
    x = cos(pi*(0:N)/N)';
    c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
    X = repmat(x,1,N+1); dX = X - X';
    D = (c*(1./c)')./(dX+eye(N+1));
    D = D - diag(sum(D,2));
end

function c = redblue_local()
    m=256; t=linspace(0,1,m)';
    c=[min(1,2*t), 1-abs(2*t-1), min(1,2*(1-t))];
end
