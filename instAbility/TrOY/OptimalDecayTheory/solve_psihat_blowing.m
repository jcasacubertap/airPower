%% solve_psihat_blowing.m
% Blowing-suction (wall-optimised) counterpart of solve_psihat_mode.m.
%
% Same interior equation as the sealed case (1.11), but the no-penetration
% wall condition is REPLACED by the natural condition that emerges when the
% wall v' is itself optimised. Optimising over the wall blowing turns the
% Dirichlet condition psihat=0 into the Neumann condition psihat'=0:
%
%        psihat'' - beta^2 psihat = i beta (duB/dy) uhat,   psihat'(+-1) = 0.
%
% At the wall this means  w' = -d psihat/dy = 0  (no spanwise slip) while
% v' = d psi/dz is free -- i.e. a pure blowing/suction wall. The OPTIMAL
% blowing profile to prescribe in a stability solver is the wall trace
%
%        v'(+-1,z) = Re[ i beta psihat(+-1) e^{i beta z} ],   w'(+-1,z) = 0.
%
% Streak mode:  u'(y,z) = Re[ uhat(y) e^{i beta z} ],  uhat = A(y) e^{i phi}.
% Base flow (default): plane Poiseuille  uB = 1 - y^2,  duB/dy = -2y.

clear; clc;

%% ----- inputs ------------------------------------------------------------
beta = 2.0;                 % spanwise wavenumber
N    = 100;                 % Chebyshev resolution in y

Afun   = @(y) (1-y.^2).*exp(-((y+0.7)/0.25).^2);     % streak amplitude  A(y)
phifun = @(y) 1.4*y;          % streak phase      phi(y)   (e.g. 0.8*y to tilt)
dUBdy  = @(y) -2*y;         % base-flow shear   duB/dy

%% ----- discretisation and assembly --------------------------------------
[D,y] = cheb(N);
D2    = D*D;
Up    = dUBdy(y);
uhat  = Afun(y).*exp(1i*phifun(y));

L = D2 - beta^2*eye(N+1);                    % operator  d^2/dy^2 - beta^2
r = 1i*beta*(Up.*uhat);                      % rhs  i beta (duB/dy) uhat

% ---- Neumann boundary conditions  psihat'(+-1) = 0  (the ONLY change) ----
L(1,:)   = D(1,:);   r(1)   = 0;             % d psihat/dy = 0 at y = +1
L(N+1,:) = D(N+1,:); r(N+1) = 0;             % d psihat/dy = 0 at y = -1

psihat = L\r;                                % solve (blowing optimum)

%% ----- cross-stream velocity (eq. 1.8) ----------------------------------
vhat =  1i*beta*psihat;                      % vhat = i beta psihat
what = -D*psihat;                            % what = -d psihat/dy  (=0 at walls)

%% ----- the OPTIMAL WALL BLOWING to prescribe in the solver --------------
vtop = 1i*beta*psihat(1);                    % vhat at y = +1
vbot = 1i*beta*psihat(N+1);                  % vhat at y = -1
vwall = @(side,z) real( ((side>0)*vtop + (side<0)*vbot).*exp(1i*beta*z) );

%% ----- verification ------------------------------------------------------
cont   = D*vhat + 1i*beta*what;              % continuity  (machine zero)
omx    = D*what - 1i*beta*vhat;              % streamwise vorticity, modal
omx_id = omx + Up.*(1i*beta*uhat);           % omega_x = -(duB/dy) du'/dz -> ~0
[~,wq] = clencurt(N); wq = wq(:);
normPg = sqrt( sum( wq.*(beta^2*abs(psihat).^2 + abs(D*psihat).^2) ) );  % ||Pi g||
fprintf('max |continuity|             : %.2e\n', max(abs(cont)));
fprintf('max |omega_x + U'' du''/dz|    : %.2e\n', max(abs(omx_id(2:N))));
fprintf('|d psihat/dy| at walls (Neu) : %.2e , %.2e\n', abs(what(1)), abs(what(N+1)));
fprintf('wall blowing vhat(+1),(-1)   : %+.4f , %+.4f\n', real(vtop), real(vbot));
fprintf('||Pi g|| (transfer at E0=1/2): %.4f\n', normPg);

%% ----- reconstruct physical fields over two wavelengths -----------------
Lz = 2*pi/beta; Nz = 240; z = linspace(0,2*Lz,Nz);
E  = exp(1i*beta*z); rec = @(a) real(a*E);
up  = rec(uhat);  psi = rec(psihat);  vp = rec(vhat);  wp = rec(what);

%% ----- plot: streak + roller streamlines + wall blowing arrows ----------
[Z,Y] = meshgrid(z,y);
figure('Color','w','Position',[100 100 820 440]);
pcolor(Z,Y,up); shading interp; hold on
colormap(redblue_local()); cmax=max(abs(up(:))); caxis([-cmax cmax]);
cb=colorbar; cb.Label.String='u'' (streak)';
contour(Z,Y,psi,12,'k','LineWidth',0.5);
iy=1:max(1,round(N/16)):N+1; iz=1:max(1,round(Nz/18)):Nz;
quiver(Z(iy,iz),Y(iy,iz),wp(iy,iz),vp(iy,iz),0.9,'k');
% mark the wall blowing along the top and bottom edges
zb = linspace(0,2*Lz,40);
quiver(zb, 1+0*zb, 0*zb, vwall(+1,zb), 0.1,'b','LineWidth',1);   % top wall v'
quiver(zb,-1+0*zb, 0*zb, vwall(-1,zb), 0.1,'b','LineWidth',1);   % bottom wall v'
xlabel('z'); ylabel('y'); axis tight
title(sprintf('Blowing-suction optimum, \\beta = %.2f (blue: wall v'')',beta));
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

function [x,w] = clencurt(N)
    theta = pi*(0:N)'/N; x = cos(theta);
    w = zeros(1,N+1); ii = 2:N; v = ones(N-1,1);
    if mod(N,2)==0
        w(1) = 1/(N^2-1); w(N+1) = w(1);
        for k = 1:N/2-1, v = v - 2*cos(2*k*theta(ii))/(4*k^2-1); end
        v = v - cos(N*theta(ii))/(N^2-1);
    else
        w(1) = 1/N^2; w(N+1) = w(1);
        for k = 1:(N-1)/2, v = v - 2*cos(2*k*theta(ii))/(4*k^2-1); end
    end
    w(ii) = 2*v/N;
end

function c = redblue_local()
    m=256; t=linspace(0,1,m)';
    c=[min(1,2*t), 1-abs(2*t-1), min(1,2*(1-t))];
end
