function DD = FD1_uneven(x, y, order)
%% Description
% Computes first derivatives on a non-uniform grid using finite differences.
% Derivative is taken along columns.
%
% Inputs:
%   x     - Grid coordinates (column vector or matrix)
%   y     - Function values at x (column vector or matrix)
%   order - Order of accuracy (2, 4, or 6). Default: 2
%
% Output:
%   DD - First derivative of y with respect to x
%
% Example:
%   dfdx = FD1d_uneven(x, f, 4);

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
        error('FD1d_uneven:invalidOrder', ...
            'Order must be 2, 4, or 6. Got %d', order);
end

if nrows < minRows
    error('FD1d_uneven:insufficientPoints', ...
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

