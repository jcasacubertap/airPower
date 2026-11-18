function DD = deriv2_uneven(x, y, order)
%% License (GNU GENERAL PUBLIC LICENSE v3)
%                  Delft Harmonic Navier-Stokes Solver
%     Copyright (C) 2023 S.H.J. Westerbeek, S. Hulshoff, H. Schuttelaars
%                          & M. Kotsonis
%     Modified (C) 2026 P.C. Boldini, Y. Chen & M. Kotsonis
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
%     Please cite this code and the paper if you have used this for your
%     publication:

% DOI:

%% Description
% Second derivative on a non-uniform grid (order 2, 4, or 6).
%
%   DD = deriv2_uneven(x, y, order)
%
%   Differentiates y twice with respect to x along columns using order-2/4/6
%   finite differences on a non-uniform grid (central in the interior, one-sided
%   at the ends). For a matrix y, each column is differentiated independently.
%
%   x holds the coordinates and y the values, either a column vector or a matrix;
%   order is 2, 4 or 6 (default 2). DD is d2y/dx2, the same size as y.

%% Inputs
% Name  size  unit   explanation
% x     (n,1) [-]    coordinate, need not be uniformly spaced
% y     (n,m) [-]    data to differentiate, column by column
% order (1,1) [-]    finite-difference order: 2, 4 or 6

%% Outputs
% Name size  unit   explanation
% DD   (n,m) [-]    second derivative of y with respect to x, one-sided at the boundaries

%% Default order
if nargin < 3
    order = 2;
end

%% Input validation
[nrows, ncols] = size(y);

switch order
    case 2
        minRows = 3;
    case 4
        minRows = 5;
    case 6
        minRows = 7;
    otherwise
        error('deriv2_uneven:invalidOrder', ...
            'Order must be 2, 4, or 6. Got %d', order);
end

if nrows < minRows
    error('deriv2_uneven:insufficientPoints', ...
        'Order %d requires at least %d rows, got %d', order, minRows, nrows);
end

%% Compute grid spacing
dx = diff(x);

%% Preallocate output
DD = zeros(size(y));

%% Compute derivatives based on order
switch order
    case 2
        % Central differences (interior nodes)
        DD(2:end-1,:) = 2 * (   dx(1:end-1,:)              .* y(3:end,:) ...
                             - (dx(1:end-1,:)+dx(2:end,:)) .* y(2:end-1,:) ...
                             +  dx(2:end,:)                .* y(1:end-2,:) ) ...
                          ./ ( dx(1:end-1,:) .* dx(2:end,:) .* (dx(1:end-1,:) + dx(2:end,:)) );
        
        % Forward differences (first node)
        DD(1,:) = 2 * (   dx(1,:)            .* y(3,:) ...
                       - (dx(1,:)+dx(2,:))   .* y(2,:) ...
                       +  dx(2,:)            .* y(1,:) ) ...
                    ./ ( dx(1,:) .* dx(2,:) .* (dx(1,:) + dx(2,:)) );
        
        % Backward differences (last node)
        DD(end,:) = 2 * (   dx(end-1,:)              .* y(end,:) ...
                         - (dx(end-1,:)+dx(end,:))   .* y(end-1,:) ...
                         +  dx(end,:)                .* y(end-2,:) ) ...
                      ./ ( dx(end-1,:) .* dx(end,:) .* (dx(end-1,:) + dx(end,:)) );
    
    case 4
        % Precompute coefficients for all interior nodes
        for i = 3:nrows-2
            % Grid positions relative to point i
            xm2 = x(i-2) - x(i);
            xm1 = x(i-1) - x(i);
            x0  = 0;
            xp1 = x(i+1) - x(i);
            xp2 = x(i+2) - x(i);
            
            % Vandermonde matrix
            A = [1,     1,     1,     1,     1;
                 xm2,   xm1,   x0,    xp1,   xp2;
                 xm2^2, xm1^2, x0^2,  xp1^2, xp2^2;
                 xm2^3, xm1^3, x0^3,  xp1^3, xp2^3;
                 xm2^4, xm1^4, x0^4,  xp1^4, xp2^4];
            
            b = [0; 0; 2; 0; 0];
            c = A \ b;
            
            DD(i,:) = c(1) * y(i-2,:) + c(2) * y(i-1,:) + c(3) * y(i,:) + c(4) * y(i+1,:) + c(5) * y(i+2,:);
        end
        
        % Node 1 (forward stencil using points 1,2,3,4,5)
        xp = x(1:5) - x(1);
        A = [ones(1,5); xp'; (xp.^2)'; (xp.^3)'; (xp.^4)'];
        c = A \ [0; 0; 2; 0; 0];
        DD(1,:) = c(1)*y(1,:) + c(2)*y(2,:) + c(3)*y(3,:) + c(4)*y(4,:) + c(5)*y(5,:);
        
        % Node 2 (forward-biased stencil using points 1,2,3,4,5)
        xp = x(1:5) - x(2);
        A = [ones(1,5); xp'; (xp.^2)'; (xp.^3)'; (xp.^4)'];
        c = A \ [0; 0; 2; 0; 0];
        DD(2,:) = c(1)*y(1,:) + c(2)*y(2,:) + c(3)*y(3,:) + c(4)*y(4,:) + c(5)*y(5,:);
        
        % Node end-1 (backward-biased stencil using points end-4,end-3,end-2,end-1,end)
        xp = x(end-4:end) - x(end-1);
        A = [ones(1,5); xp'; (xp.^2)'; (xp.^3)'; (xp.^4)'];
        c = A \ [0; 0; 2; 0; 0];
        DD(end-1,:) = c(1)*y(end-4,:) + c(2)*y(end-3,:) + c(3)*y(end-2,:) + c(4)*y(end-1,:) + c(5)*y(end,:);
        
        % Node end (backward stencil using points end-4,end-3,end-2,end-1,end)
        xp = x(end-4:end) - x(end);
        A = [ones(1,5); xp'; (xp.^2)'; (xp.^3)'; (xp.^4)'];
        c = A \ [0; 0; 2; 0; 0];
        DD(end,:) = c(1)*y(end-4,:) + c(2)*y(end-3,:) + c(3)*y(end-2,:) + c(4)*y(end-1,:) + c(5)*y(end,:);

    case 6
        % 7-point stencils. Weights solve a Vandermonde system scaled by h0 so the
        % matrix entries stay O(1) (well-conditioned for any spacing): in u = xp/h0,
        % A*ctil = [0;0;2;0;0;0;0], then c = ctil/h0^2 for the second derivative.
        np = 7;  bvec = [0;0;2;0;0;0;0];
        % Interior nodes: centered 7-point stencil
        for i = 4:nrows-3
            xp = x(i-3:i+3) - x(i);  h0 = max(abs(xp));
            c  = (((xp(:).'/h0).^((0:np-1)')) \ bvec) / h0^2;
            DD(i,:) = c.' * y(i-3:i+3,:);
        end
        % Leading boundary nodes: forward-biased 7-point stencil (points 1..7)
        for i = 1:3
            xp = x(1:np) - x(i);  h0 = max(abs(xp));
            c  = (((xp(:).'/h0).^((0:np-1)')) \ bvec) / h0^2;
            DD(i,:) = c.' * y(1:np,:);
        end
        % Trailing boundary nodes: backward-biased 7-point stencil (points end-6..end)
        for i = nrows-2:nrows
            xp = x(nrows-np+1:nrows) - x(i);  h0 = max(abs(xp));
            c  = (((xp(:).'/h0).^((0:np-1)')) \ bvec) / h0^2;
            DD(i,:) = c.' * y(nrows-np+1:nrows,:);
        end
end

end

