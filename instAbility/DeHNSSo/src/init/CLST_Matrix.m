function [Lt,Lx,Ly,Lz,Lq,Vxx,Vxy,Vyy,Vxz,Vyz,Vzz] = CLST_Matrix(Re, Pr, Ec, ...
    rk, ryk, tk, tyk, tyyk, uk, uyk, uyyk, wk, wyk, ...
    muk, mutk, muttk, muyk, pk, eT, gM2)
%CLST_MATRIX  Assemble local stability matrices for one grid point.
%
%   [Lt,Lx,Ly,Lz,Lq,Vxx,Vxy,Vyy,Vxz,Vyz,Vzz] = CLST_Matrix(...)
%
%   Returns the 5x5 coefficient matrices for the compressible LST
%   eigenvalue problem with variables (rho, u, v, w, T).
%   Ideal gas EOS: p = rho*T/(gamma*Ma^2), e = cv*T.
%
%   Inputs:
%     Re, Pr, Ec  — Reynolds, Prandtl, Eckert numbers
%     rk, ryk     — base density and its y-derivative
%     tk, tyk, tyyk — base temperature and derivatives
%     uk, uyk, uyyk — base streamwise velocity and derivatives
%     wk, wyk     — base spanwise velocity and derivative
%     muk, mutk, muttk — viscosity and temperature derivatives
%     muyk        — y-derivative of viscosity (= mutk * tyk)
%     pk          — base pressure = rk*tk/gM2
%     eT          — de/dT = cv = 1/(gamma*(gamma-1)*Ma^2)
%     gM2         — gamma*Ma^2

% EOS derivatives (ideal gas: p = rho*T/gM2)
p_rho = tk / gM2;
p_T   = rk / gM2;
p_rhoT = 1 / gM2;

%% Lt (temporal: -iw coefficient)
Lt = diag([1, rk, rk, rk, rk*eT]);

%% Lx (streamwise: ia coefficient)
Lx = zeros(5);
Lx(1,1) = uk;              Lx(1,2) = rk;
Lx(2,1) = p_rho;           Lx(2,2) = rk*uk;
Lx(2,3) = -muyk/Re;        Lx(2,5) = p_T;
Lx(3,2) = 2/3*muyk/Re;     Lx(3,3) = rk*uk;
Lx(3,5) = -mutk*uyk/Re;
Lx(4,4) = rk*uk;
Lx(5,2) = pk;              Lx(5,3) = -2*muk*uyk/Re;
Lx(5,5) = rk*uk*eT;

%% Ly (wall-normal: d/dy coefficient)
Ly = zeros(5);
Ly(1,3) = rk;
Ly(2,2) = -muyk/Re;        Ly(2,5) = -mutk*uyk/Re;
Ly(3,1) = p_rho;            Ly(3,3) = -4/3*muyk/Re;
Ly(3,5) = p_T;
Ly(4,4) = -muyk/Re;
Ly(5,2) = -2*muk*uyk/Re;   Ly(5,3) = pk;
Ly(5,5) = -(muyk + mutk*tyk) / (Re*Pr*Ec);

%% Lz (spanwise: ib coefficient)
Lz = zeros(5);
Lz(1,4) = rk;
Lz(3,4) = 2/3*muyk/Re;
Lz(4,1) = p_rho;           Lz(4,3) = -muyk/Re;
Lz(4,5) = p_T;
Lz(5,4) = pk;

%% Lq (zero-derivative terms)
Lq = zeros(5);
Lq(1,3) = ryk;
Lq(2,3) = rk*uyk;
Lq(2,5) = -(mutk*uyyk + muttk*uyk*tyk) / Re;
Lq(3,1) = p_rhoT*tyk;      Lq(3,5) = p_rhoT*ryk;
Lq(4,3) = rk*wyk;
Lq(5,3) = rk*eT*tyk;
Lq(5,5) = -(tyyk*mutk + tyk^2*muttk)/(Re*Pr*Ec) ...
          - mutk*(uyk^2 + wyk^2)/Re;

%% Viscous second-derivative matrices
Vxx = diag([0, -4/3*muk/Re, -muk/Re,     -muk/Re,     -muk/(Re*Pr*Ec)]);
Vyy = diag([0, -muk/Re,     -4/3*muk/Re,  -muk/Re,     -muk/(Re*Pr*Ec)]);
Vzz = diag([0, -muk/Re,     -muk/Re,      -4/3*muk/Re, -muk/(Re*Pr*Ec)]);

Vxy = zeros(5);  Vxy(2,3) = -1/3*muk/Re;  Vxy(3,2) = -1/3*muk/Re;
Vxz = zeros(5);  Vxz(2,4) = -1/3*muk/Re;  Vxz(4,2) = -1/3*muk/Re;
Vyz = zeros(5);  Vyz(3,4) = -1/3*muk/Re;  Vyz(4,3) = -1/3*muk/Re;

end
