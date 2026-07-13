function DD = FD2d_uneven(x, y, order)
%% Description
% Computes second derivatives on a non-uniform grid using finite differences.
% Derivative is taken along columns.
%
% Inputs:
%   x     - Grid coordinates (column vector or matrix)
%   y     - Function values at x (column vector or matrix)
%   order - Order of accuracy (2 or 4). Default: 2
%
% Output:
%   DD - Second derivative of y with respect to x
%
% Example:
%   d2fdx2 = FD2d_uneven(x, f, 4);

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
    otherwise
        error('FD2d_uneven:invalidOrder', ...
            'Order must be 2 or 4. Got %d', order);
end

if nrows < minRows
    error('FD2d_uneven:insufficientPoints', ...
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
end

end

