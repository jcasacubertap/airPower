function [EV,u,v,w,p]=ILST_solver(Re,Ur,Wr,dyUr,dyWr,omega,beta,N_LST,dzdy,d2zdy2,bc_top,alpha_shift)
%ILST_SOLVER  Solve the local stability eigenvalue problem (Orr-Sommerfeld/Squire).
%
%   [EV, u, v, w, p] = ILST_solver(Re, Ur, Wr, dyUr, dyWr, omega, beta, N_LST, dzdy, d2zdy2)
%
%   Inputs:
%     Re      — Reynolds number at the inflow station
%     Ur, Wr  — base-flow velocity profiles [N_LST x 1]
%     dyUr, dyWr — wall-normal derivatives of base flow [N_LST x 1]
%     omega   — angular frequency
%     beta    — spanwise wavenumber
%     N_LST   — number of Chebyshev collocation points
%     dzdy, d2zdy2 — Malik mapping derivatives [N_LST x 1]
%
%   Outputs:
%     EV      — eigenvalue spectrum (streamwise wavenumbers alpha)
%     u,v,w,p — eigenvectors [N_LST x n_eig] for each eigenvalue

    % Chebyshev operators
    Dn  = chebDiff(N_LST);
    Dn2 = Dn * Dn;
    II  = eye(4);
    Total_Dn  = kron(Dn, II);
    Total_DDn = kron(Dn2, II);

    % Pre-allocate block-diagonal matrices
    n4 = 4 * N_LST;
    AAL1 = zeros(n4);
    AAL2 = zeros(n4);
    AAL3 = zeros(n4);
    BBR1 = zeros(n4);
    BBR2 = zeros(n4);

    % Incompressible constants (mu_y=0 simplifies many terms)
    mu = 1;  Rh = 1;  factp = 1;

    % Viscous matrices (constant, diagonal — compute once)
    Vxx = diag([0, -4/3*mu/Re, -mu/Re, -mu/Re]);
    Vyy = diag([0, -mu/Re, -4/3*mu/Re, -mu/Re]);
    Vzz = diag([0, -mu/Re, -mu/Re, -4/3*mu/Re]);
    Vxy = zeros(4); Vxy(2,3) = -1/3*mu/Re; Vxy(3,2) = -1/3*mu/Re;
    Vxz = zeros(4); Vxz(2,4) = -1/3*mu/Re; Vxz(4,2) = -1/3*mu/Re;
    Vyz = zeros(4); Vyz(3,4) = -1/3*mu/Re; Vyz(4,3) = -1/3*mu/Re;

    % Build block-diagonal matrices (vectorised where possible)
    for i = 1:N_LST
        % Variable matrices (depend on U, W, U_y, W_y at this point)
        Lt = diag([0, 1, 1, 1]);

        Lx = zeros(4);
        Lx(1,2) = Rh;  Lx(2,1) = 1/factp;  Lx(2,2) = Rh*Ur(i);
        Lx(3,3) = Rh*Ur(i);  Lx(4,4) = Rh*Ur(i);

        Ly = zeros(4);
        Ly(1,3) = Rh;  Ly(3,1) = 1/factp;

        Lz = zeros(4);
        Lz(1,4) = Rh;  Lz(2,2) = Rh*Wr(i);  Lz(3,3) = Rh*Wr(i);
        Lz(4,1) = 1/factp;  Lz(4,4) = Rh*Wr(i);

        Lq = zeros(4);
        Lq(2,3) = Rh*dyUr(i);  Lq(4,3) = Rh*dyWr(i);

        % Combine
        A1 = -1i*omega*Lt + 1i*beta*Lz + Lq - beta^2*Vzz;
        A2 = Ly + 1i*beta*Vyz;
        A3 = Vyy;
        B1 = -1i*Lx + beta*Vxz;
        B2 = -1i*Vxy;

        % Apply Malik mapping
        dz = dzdy(i);
        d2z = d2zdy2(i);

        idx = (i-1)*4+1 : i*4;
        AAL1(idx,idx) = A1;
        AAL2(idx,idx) = A2*dz + A3*d2z;
        AAL3(idx,idx) = A3*dz^2;
        BBR1(idx,idx) = B1;
        BBR2(idx,idx) = B2*dz;

        %Antonis 2026% 
        % idx = (i-1)*4+1 : i*4;
        % AAL1(idx,idx) = A1;
        % AAL2(idx,idx) = A2;
        % AAL3(idx,idx) = A3;
        % BBR1(idx,idx) = B1;
        % BBR2(idx,idx) = B2;


    end

    AA = AAL1 + AAL2 * Total_Dn + AAL3 * Total_DDn;
    BB = BBR1 + BBR2 * Total_Dn;

    % Wall BCs — always Dirichlet: u=v=w=0 at y=0 (no-slip)
    % Inputs are flipped so index 1 = wall, index N_LST = freestream.
    AA(2:4,:) = 0;
    AA(2,2) = 1;  AA(3,3) = 1;  AA(4,4) = 1;
    BB(2:4,:) = 0;

    % Farfield BCs — H_DR (u=v=w=0) or H_NM (d/dy=0) at y=H
    if ~exist('bc_top','var'); bc_top = 'H_DR'; end
    AA(n4-2:n4,:) = 0;
    BB(n4-2:n4,:) = 0;
    if strcmpi(bc_top,'H_NM')
        % d/dz = 0 at freestream: use the corresponding rows of Total_Dn
        AA(n4-2,:) = Total_Dn(n4-2,:);
        AA(n4-1,:) = Total_Dn(n4-1,:);
        AA(n4,  :) = Total_Dn(n4,  :);
    else
        % Dirichlet: u=v=w=0 at freestream
        AA(n4-2,n4-2) = 1;  AA(n4-1,n4-1) = 1;  AA(n4,n4) = 1;
    end

    % Solve generalised eigenvalue problem
    if exist('alpha_shift','var') && ~isempty(alpha_shift)
        % Mode continuation: eigs with shift (fast, one EV)
        opts_eigs.tol   = 1e-10;
        opts_eigs.maxit = 300;
        opts_eigs.v0    = ones(n4, 1);
        [Evec_s, EV_s] = eigs(AA, BB, 1, alpha_shift, opts_eigs);
        EV   = EV_s;
        Evec = Evec_s;
    else
        [Evec, P] = eig(AA, BB);
        EV = diag(P);
    end

    % Extract eigenvectors
    p = flip(Evec(1:4:end,:));
    u = flip(Evec(2:4:end,:));
    v = flip(Evec(3:4:end,:));
    w = flip(Evec(4:4:end,:));
end
