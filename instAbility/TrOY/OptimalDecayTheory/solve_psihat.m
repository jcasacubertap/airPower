%% solve_psihat.m
% Optimal de-streaking cross-stream field for a single spanwise Fourier mode.
% Solves the boundary-value problem for the cross-plane streamfunction,
%
%        d^2 psihat/dy^2 - beta^2 psihat = i beta (duB/dy) uhat,
%
% for a prescribed streak mode  u'(y,z) = Re[ uhat(y) e^{i beta z} ],
% uhat(y) = A(y) e^{i phi(y)}.  The cross-stream velocity follows as
% vhat = i beta psihat,  what = -d psihat/dy.
%
% The two cases differ ONLY in the wall condition:
%
%   wallBC = 'sealed'   ->  psihat(+-1) = 0          (Dirichlet)
%                           no penetration: v' = 0 at the wall, w' free.
%                           The optimal interior roller.
%
%   wallBC = 'blowing'  ->  d psihat/dy (+-1) = 0    (Neumann)
%                           the natural condition when the wall v' is itself
%                           optimised: w' = 0 at the wall (no spanwise slip),
%                           v' free. The wall trace
%                              v'(+-1,z) = Re[ i beta psihat(+-1) e^{i beta z} ]
%                           is the OPTIMAL blowing/suction to prescribe in a
%                           stability solver.
%
% Discretisation: Chebyshev collocation in y on [-1,1]. A single wavenumber
% beta is solved; a general spanwise-periodic streak is recovered by calling
% this for each constituent beta and superposing.
%
% Base flow (default): plane Poiseuille  uB = 1 - y^2,  duB/dy = -2y.

clear; clc;

%% ----- inputs ------------------------------------------------------------
wallBC = 'blowing';         % 'sealed' | 'blowing'

beta = 2.0;                 % spanwise wavenumber
N    = 100;                 % Chebyshev resolution in y

% Library of streak amplitudes A(y):
%   @(y) 1 - y.^2;                                 canonical, symmetric
%   @(y) sin(pi*y);                                antisymmetric
%   @(y) (1-y.^2).*exp(-((y+0.7)/0.25).^2);        near-wall packet
%   @(y) (1-y.^2).*cos(2*pi*y);                    multi-lobe
%   @(y) cos(pi*y/2);                              smooth single peak
Afun   = @(y) 1 - y.^2;     % streak amplitude  A(y)
phifun = @(y) 0.8*y;        % streak phase      phi(y)   (0 = untilted)
dUBdy  = @(y) -2*y;         % base-flow shear   duB/dy

%% ----- discretisation and assembly --------------------------------------
[D,y] = cheb(N);            % nodes y (col, +1..-1) and 1st-derivative matrix
D2    = D*D;
Up    = dUBdy(y);
uhat  = Afun(y).*exp(1i*phifun(y));          % uhat(y) = A e^{i phi}

L = D2 - beta^2*eye(N+1);                    % operator  d^2/dy^2 - beta^2
r = 1i*beta*(Up.*uhat);                      % rhs  i beta (duB/dy) uhat

%% ----- wall conditions: the ONLY difference between the two cases -------
switch lower(wallBC)
    case 'sealed'                            % psihat(+-1) = 0
        L(1,:)   = 0; L(1,1)     = 1; r(1)   = 0;    % y = +1
        L(N+1,:) = 0; L(N+1,N+1) = 1; r(N+1) = 0;    % y = -1
    case 'blowing'                           % d psihat/dy (+-1) = 0
        L(1,:)   = D(1,:);   r(1)   = 0;             % y = +1
        L(N+1,:) = D(N+1,:); r(N+1) = 0;             % y = -1
    otherwise
        error('wallBC must be ''sealed'' or ''blowing''.');
end

psihat = L\r;                                % solve

%% ----- cross-stream velocity --------------------------------------------
vhat =  1i*beta*psihat;                      % vhat = i beta psihat
what = -D*psihat;                            % what = -d psihat/dy

vtop = vhat(1);                              % wall trace at y = +1
vbot = vhat(N+1);                            % wall trace at y = -1
vwall = @(side,z) real( ((side>0)*vtop + (side<0)*vbot).*exp(1i*beta*z) );

%% ----- verification ------------------------------------------------------
cont   = D*vhat + 1i*beta*what;              % continuity  (machine zero)
omx    = D*what - 1i*beta*vhat;              % streamwise vorticity, modal
omx_id = omx + Up.*(1i*beta*uhat);           % omega_x = -(duB/dy) du'/dz -> ~0
[~,wq] = clencurt(N); wq = wq(:);
normPg = sqrt( sum( wq.*(beta^2*abs(psihat).^2 + abs(D*psihat).^2) ) );

fprintf('--- wallBC = %s,  beta = %.2f ---\n', wallBC, beta);
fprintf('max |continuity|             : %.2e\n', max(abs(cont)));
fprintf('max |omega_x + U'' du''/dz|    : %.2e\n', max(abs(omx_id(2:N))));
fprintf('|psihat|    at walls         : %.2e , %.2e\n', abs(psihat(1)), abs(psihat(N+1)));
fprintf('|d psihat/dy| at walls       : %.2e , %.2e\n', abs(what(1)), abs(what(N+1)));
fprintf('wall v-hat  (+1) , (-1)      : %+.4f , %+.4f\n', real(vtop), real(vbot));
fprintf('||Pi g||  (transfer, E0=1/2) : %.4f\n', normPg);

%% ----- reconstruct physical fields over two wavelengths -----------------
Lz = 2*pi/beta; Nz = 240; z = linspace(0,2*Lz,Nz);
E  = exp(1i*beta*z); rec = @(a) real(a*E);
up  = rec(uhat);  psi = rec(psihat);  vp = rec(vhat);  wp = rec(what);

%% ----- plot: streak + roller streamlines + (w',v') quiver ---------------
[Z,Y] = meshgrid(z,y);
figure('Color','w','Position',[100 100 820 440]);
pcolor(Z,Y,up); shading interp; hold on
colormap(redblue_local()); cmax = max(abs(up(:))); caxis([-cmax cmax]);
cb = colorbar; cb.Label.String = 'u'' (streak)';
contour(Z,Y,psi,12,'k','LineWidth',0.5);
iy = 1:max(1,round(N/16)):N+1; iz = 1:max(1,round(Nz/18)):Nz;
quiver(Z(iy,iz),Y(iy,iz),wp(iy,iz),vp(iy,iz),0.9,'k');

if strcmpi(wallBC,'blowing')                 % draw the wall actuation
    zb = linspace(0,2*Lz,40);
    quiver(zb, 1+0*zb, 0*zb, vwall(+1,zb), 0.1,'b','LineWidth',1);
    quiver(zb,-1+0*zb, 0*zb, vwall(-1,zb), 0.1,'b','LineWidth',1);
    ttl = sprintf('Blowing-suction optimum, \\beta = %.2f (blue: wall v'')',beta);
else
    ttl = sprintf('Sealed-wall optimum, \\beta = %.2f',beta);
end
xlabel('z'); ylabel('y'); axis tight; title(ttl); hold off

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
