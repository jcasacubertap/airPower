function DD = deriv1_uneven(x, y, order)
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
% First derivative on a non-uniform grid (order 2, 4, or 6).
%
%   DD = deriv1_uneven(x, y, order)
%
%   Differentiates y with respect to x along columns using order-2/4/6 finite
%   differences on a non-uniform grid (central in the interior, one-sided at the
%   ends). For a matrix y, each column is differentiated independently.
%
%   x holds the coordinates and y the values, either a column vector or a matrix;
%   order is 2, 4 or 6 (default 2). DD is dy/dx, the same size as y.

%% Inputs
% Name  size  unit   explanation
% x     (n,1) [-]    coordinate, need not be uniformly spaced
% y     (n,m) [-]    data to differentiate, column by column
% order (1,1) [-]    finite-difference order: 2, 4 or 6

%% Outputs
% Name size  unit   explanation
% DD   (n,m) [-]    first derivative of y with respect to x, one-sided at the boundaries

%% Default order
if nargin < 3
    order = 2;
end

%% Input validation
[nrows, ~] = size(y);

switch order
    case 2
        minRows = 3;
        nStencil = 3;
    case 4
        minRows = 5;
        nStencil = 5;
    case 6
        minRows = 7;
        nStencil = 7;
    otherwise
        error('deriv1_uneven:invalidOrder', ...
            'Order must be 2, 4, or 6. Got %d', order);
end

if nrows < minRows
    error('deriv1_uneven:insufficientPoints', ...
        'Order %d requires at least %d rows, got %d', order, minRows, nrows);
end

%% Compute grid spacing
dx = diff(x);

%% Preallocate output
DD = zeros(size(y));

%% Compute derivatives based on order
switch order
    case 2
        % Central differences
        DD(2:end-1,:) = (   dx(1:end-1,:).^2                     .* y(3:end  ,:) ...
                         + (dx(2:end  ,:).^2 - dx(1:end-1,:).^2) .* y(2:end-1,:) ...
                         -  dx(2:end  ,:).^2                     .* y(1:end-2,:) ) ...
                       ./ ( dx(1:end-1,:) .* dx(2:end,:) .* (dx(1:end-1,:) + dx(2:end,:)) );
        
        % Forward differences (first node)
        DD(1,:) = ( -dx(1,:).^2                         .* y(3,:) ...
                   + (dx(1,:) + dx(2,:)).^2             .* y(2,:) ...
                   - (dx(2,:).^2 + 2.*dx(1,:).*dx(2,:)) .* y(1,:) ) ...
                 ./ ( dx(1,:) .* dx(2,:) .* (dx(1,:) + dx(2,:)) );
        
        % Backward differences (last node)
        DD(end,:) = ( (dx(end-1,:).^2 + 2.*dx(end-1,:).*dx(end,:)) .* y(end,:) ...
                     - (dx(end-1,:) + dx(end,:)).^2                .* y(end-1,:) ...
                     +  dx(end,:).^2                               .* y(end-2,:) ) ...
                   ./ ( dx(end-1,:) .* dx(end,:) .* (dx(end-1,:) + dx(end,:)) );
    
    case 4
        nBoundary = 2;  % Number of boundary nodes on each side
        
        % Interior nodes - central 5-point stencil
        for i = (nBoundary+1):(nrows-nBoundary)
            xp = x(i-2:i+2) - x(i);
            A = [ones(1,5); xp'; (xp.^2)'; (xp.^3)'; (xp.^4)'];
            b = [0; 1; 0; 0; 0];
            c = A \ b;
            DD(i,:) = c(1)*y(i-2,:) + c(2)*y(i-1,:) + c(3)*y(i,:) + c(4)*y(i+1,:) + c(5)*y(i+2,:);
        end
        
        % Left boundary nodes (1, 2)
        for i = 1:nBoundary
            xp = x(1:5) - x(i);
            A = [ones(1,5); xp'; (xp.^2)'; (xp.^3)'; (xp.^4)'];
            b = [0; 1; 0; 0; 0];
            c = A \ b;
            DD(i,:) = c(1)*y(1,:) + c(2)*y(2,:) + c(3)*y(3,:) + c(4)*y(4,:) + c(5)*y(5,:);
        end
        
        % Right boundary nodes (end-1, end)
        for i = (nrows-nBoundary+1):nrows
            xp = x(end-4:end) - x(i);
            A = [ones(1,5); xp'; (xp.^2)'; (xp.^3)'; (xp.^4)'];
            b = [0; 1; 0; 0; 0];
            c = A \ b;
            DD(i,:) = c(1)*y(end-4,:) + c(2)*y(end-3,:) + c(3)*y(end-2,:) + c(4)*y(end-1,:) + c(5)*y(end,:);
        end
    
    case 6
        nBoundary = 3;  % Number of boundary nodes on each side
        
        % Interior nodes - central 7-point stencil
        for i = (nBoundary+1):(nrows-nBoundary)
            xp = x(i-3:i+3) - x(i);
            A = [ones(1,7); xp'; (xp.^2)'; (xp.^3)'; (xp.^4)'; (xp.^5)'; (xp.^6)'];
            b = [0; 1; 0; 0; 0; 0; 0];
            c = A \ b;
            DD(i,:) = c(1)*y(i-3,:) + c(2)*y(i-2,:) + c(3)*y(i-1,:) + c(4)*y(i,:) + ...
                      c(5)*y(i+1,:) + c(6)*y(i+2,:) + c(7)*y(i+3,:);
        end
        
        % Left boundary nodes (1, 2, 3)
        for i = 1:nBoundary
            xp = x(1:7) - x(i);
            A = [ones(1,7); xp'; (xp.^2)'; (xp.^3)'; (xp.^4)'; (xp.^5)'; (xp.^6)'];
            b = [0; 1; 0; 0; 0; 0; 0];
            c = A \ b;
            DD(i,:) = c(1)*y(1,:) + c(2)*y(2,:) + c(3)*y(3,:) + c(4)*y(4,:) + ...
                      c(5)*y(5,:) + c(6)*y(6,:) + c(7)*y(7,:);
        end
        
        % Right boundary nodes (end-2, end-1, end)
        for i = (nrows-nBoundary+1):nrows
            xp = x(end-6:end) - x(i);
            A = [ones(1,7); xp'; (xp.^2)'; (xp.^3)'; (xp.^4)'; (xp.^5)'; (xp.^6)'];
            b = [0; 1; 0; 0; 0; 0; 0];
            c = A \ b;
            DD(i,:) = c(1)*y(end-6,:) + c(2)*y(end-5,:) + c(3)*y(end-4,:) + c(4)*y(end-3,:) + ...
                      c(5)*y(end-2,:) + c(6)*y(end-1,:) + c(7)*y(end,:);
        end
end

end

