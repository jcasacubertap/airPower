function [ u,v,w,y,Du,DDu,Dv,DDv,Dw,DDw,beta ] = solver_IBL_spectral( nx,ny,S,L,y_i,H,nu,Ue,We )

% This function calculates the non-similar boundary layer on the domain 
% (x,y) in [0,L]x[0,H] using nx x ny uniformly-cosine spaced gridpoints, 
% the point (x,y) = (0,0) being the leading edge, given the edge velocities
% and the kinematic viscosity. The discretization used is the implicit 
% scheme in x and Chebyshev in y, second and ny^th order accurate in x and 
% y, the backslash solver is used to solve the problem. Inflow b.c. are
% produced by a Falkner Skan Cooke solver at equivalent Hartree parameter
% to the local external velocity. The nonlinear momentum equation is
% linearised by simple iteration in the same iteration for the continuity
% equation. A combine threshold of 1e-10 for normalised u and v is used.
% First and second order streamwise derivatives are used to predict the next station,
% greatly accelerating convergence. A check for laminar separation is
% performed at each iteration. If separation is detected the code returns.
%
% --- Input ---
% S   : (streamwise) domain start
% L   : (streamwise) domain end
% y_i : (wall-normal) Chebyshev node coordinate median 
% H   : (wall-normal) domain height
% nu  : kinematic viscosity
% nx  : number of gridpoints in the streamwise direction
% ny  : number of gridpoints in the wall-normal direction
% Ue  : streamwise freestream velocity distribution (1 x nx array, increa-
%       sing index -> increasing coordinate)
% We  : spanwise freestream velocity (constant, scalar)
%
% --- Output ---
% u   : streamwise velocity component (ny x nx matrix, increasing index -> 
%       increasing coordinate)
% v   : wall-normal velocity component
% w   : spanwise velocity component
% y   : y-coordinates of the Chebyshev nodes
% D*  : Derivatives of the velocity components, only computed if queried

% the nodes' x-coordinates are obtained by: 
% x = S:(L-S)/(nx-1):L;

addpath(genpath('../../Cases/crossFlowTunnelOF/DataOutputOF/'))

fprintf('Starting IBL solver.\n');


%% convergence error threshold for both u and v

thres=1e-10;

%% Initialize fields

u = zeros(ny,nx);
v = zeros(ny,nx);
w = zeros(ny,nx);

% Check finiteness Ue at the inflow
if Ue(1) == 0
    Ue(1) = 10^-4;          % Incoming (a non-zero value is needed)
end

% Check constancy We
if ~isscalar(We)
    error(['The input freestream spanwise velocity We is not a scalar, but: [',num2str(size(We)),']']);
end

% Determine uniform grid spacing
dx = (L-S)/(nx-1); % N nodes -> N-1 line segments

% Initialize differentiation matrices
[yvec,DM] = chebdif(ny,2); % differentiation matrices
[y,D1,D2] = MappingMalik(H,y_i,yvec',DM(:,:,1),DM(:,:,2));
clear DM yvec

% Initialize identity matrix
I = eye(size(D1)); 

%% Initialize inflow with a similar solution at equivalent Hartree parameter

% Extract Falkner-Skan parameters from external pressure gradient
dUe = -FD1d4o(Ue,dx);   % high order gradient of external velocity distribution

m = (dUe(1))*S/Ue(1);       % Falkner-Skan parameter
m = min(m,1);               % Limit m to m<=1 (stagnation)

beta = 2*m/(m+1);           % Hartree parameter
beta = max(beta,-0.19884);  % Limit beta to beta>=-0.19884 (separation)

%fprintf(' Using Falkner-Skan inflow profile with beta=%.3f.\n',beta)
eta_inf = 30; % Similar flow freestream condition is invariant to the case
eta_n = 2000;

% Solve FSC in real coordinates
fprintf('Using Falkner-Skan inflow profile with beta=%.3f.\n',beta)
[y_init,u_init,v_init,w_init] = FalknerSkanCookeReal(eta_inf,eta_n,beta,nu,S,Ue(1),We,true);

% Interpolate solution on collocation points grid
u(:,1) = interp1(y_init,u_init,y,'linear',u_init(end))';
v(:,1) = interp1(y_init,v_init,y,'linear','extrap')';
w(:,1) = interp1(y_init,w_init,y,'linear',w_init(end))';

% Duplicate profile at second station
u(:,2) = u(:,1);
v(:,2) = v(:,1);
w(:,2) = w(:,1);

%% Initialize LHS & RHS including BCs
LHS (1,1)   = 1; % Dirichlet type BC at freestream
LHS (ny,ny) = 1; % Dirichlet type BC at wall
RHSu(ny,1)  = 0; % No-slip u 
RHSw(ny,1)  = 0; % No-slip w

%% Start up Inflow Solution (1st order accurate in x)

tic
i = 1;     % profile x-index
in = 2:ny-1; 

errv=1;
err=1;
count=0;
while errv>thres && err>thres     % continuity iteration

    count=count+1;

    %% Momentum
    
    % Build LHS
    LHS(in,:) = diag(u(in,i+1)/dx)*I(in,:) + diag(v(in,i+1))*D1(in,:) - nu*D2(in,:);
    % Build RHS u component
    RHSu(in)  =      u(in,i+1).*u(in,i)/dx + Ue(i+1)*(Ue(i+1)-Ue(i))/dx;
    % Freestream conditions
    RHSu(1)  = Ue(i+1); 
    % Solve problems
    unew     = LHS\RHSu;
    % Convergence error 
    err = max(abs( unew-u(:,i+1) )) / Ue(i+1);
    % Update u
    u(:,i+1) = unew;
    
    %% Continuity
    
    vold=v(:,i+1);
    % Integrate v
    v(:,i+1) = cumtrapz(y,-(u(:,i+1)-u(:,i))/dx);
    v(:,i+1) = (v(:,i+1)-v(end,i+1));
    % Convergence error
    errv = max(abs( vold-v(:,i+1) )) / Ue(i+1);
    
end
fprintf(' Iterations for IC: %u. Remaining error: %3.2g\n',count,err);

% Check for laminar separation
du = (D1*u(:,i+1))/(Ue(i+1)/H);

if du(end)<0
    fprintf('\n Separating inflow detected!\n Ending BL analysis.\n')
    u = u(:,1:2); 
    v = v(:,1:2);

    % Compute w component
    RHSw(in) = u(in,i+1).*w(in,i)/dx;
    RHSw(1) = We; 
    w(:,i+1) = LHS\RHSw;
    w = w(:,1:2);
    return
end

% Compute w component after convergence
RHSw(in) = u(in,i+1).*w(in,i)/dx;
RHSw(1) = We; 
w(:,i+1) = LHS\RHSw;

% Estimate next station using finite difference derivatives (first order)
v(:,3) = 2*v(:,2)-v(:,1);
u(:,3) = 2*u(:,2)-u(:,1);


%% March in the x direction  (second order in x)
breaker = false;
percentage = 1;
fprintf(' Starting marching method.')
fprintf(' Progress: %.0f percent',0)
for i = 2:nx-1   
    
    errv=1;     % start-up value
    err=1;
    count=0;   % iteration counter

    while errv>thres && err>thres 
        count=count+1;

        %% Momentum
        
        % Build LHS
        LHS(in,:) = diag(3*u(in,i+1)/(2*dx))*I(in,:) + diag(v(in,i+1))*D1(in,:) - nu*D2(in,:);        
        % Build RHS u component
        RHSu(in)  =        u(in,i+1).*(4*u(in,i)-u(in,i-1))/(2*dx) ...
                        + Ue(i+1)*(3*Ue(i+1)-4*Ue(i)+Ue(i-1))/(2*dx);
        % Freestream conditions
        RHSu(1) = Ue(i+1);
        % Solve problems
        unew     = LHS\RHSu;
        % Convergence error
        err = max(abs( unew-u(:,i+1) )) / Ue(i+1);
        % Update u and error
        u(:,i+1) = unew; 
        
        %% Continuity
        
        vold=v(:,i+1);
        % Integrate v
        v(:,i+1) = cumtrapz(y,-(3*u(:,i+1)-4*u(:,i)+u(:,i-1))/(2*dx));
        v(:,i+1) = v(:,i+1)-v(end,i+1);
        % Convergence error
        errv = max(abs( v(:,i+1)-vold )) / Ue(i+1);
    
        %% Check for separation

        du = (D1*u(:,i+1))/(Ue(i+1)/H); % first derivative
    
        if du(end)<0
            fprintf('\n Separation detected at x/c=%.2f with Ue=%.3f.\n Ending BL analysis.\n',i/nx,Ue(i))
            % Toggle stop of marching (outer loop)
            breaker = true;
            break
        end
    end
    
    % Compute w component after convergence
    RHSw(in) = u(in,i+1).*(4*w(in,i)-w(in,i-1))/(2*dx);
    RHSw(1) = We; 
    w(:,i+1) = LHS\RHSw;
    
    % Stop marching if separation was found
    if breaker
        % Select relevant part of the solution
        nx=i;
        u = u(:,1:i); 
        v = v(:,1:i);
        w = w(:,1:i);
        break
    end
    
    % Notify progress
    if 100*(i+1)/nx == percentage
        for j = 1:ceil(log10(percentage))+9; fprintf('\b'); end
        fprintf('%.0f percent ',percentage)
        percentage = percentage + 1;
    end
    
    % Estimate next station usinge finite difference derivatives (second order)
    v(:,i+2) = (5*v(:,i+1)-4*v(:,i)+v(:,i-1))/2;
    u(:,i+2) = (5*u(:,i+1)-4*u(:,i)+u(:,i-1))/2;
    
end

% Remove next station estimation if full solution was obtained
if i==nx-1
    v=v(:,1:end-1);
    u=u(:,1:end-1);
    fprintf('\n')
end

fprintf(' Computation concluded in %.2fs!\n',toc)

%% Outputs

% Produce additional derivatives only if requested
if nargout > 4
    for i=1:nx
        Du(:,i)=D1*u(:,i);
        DDu(:,i)=D2*u(:,i);
        Dv(:,i)=D1*v(:,i);
        DDv(:,i)=D2*v(:,i);
        Dw(:,i)=D1*w(:,i);
        DDw(:,i)=D2*w(:,i);    
    end
end

end

