function Grid = compute_cheb_operators_channel(Grid)
%COMPUTE_CHEB_OPERATORS  Chebyshev D1, D2 with Malik wall-normal mapping.
%
%   Malik mapping:  z = -1 + (1+b)*eta / (a + eta)
%   where  a = H*y_i / (H - 2*y_i),  b = 1 + 2*a/H
%
%   D1 = diag(dz/deta) * Dn        ->  d/d(eta)
%   D2 = diag(d2z/deta2)*Dn + diag((dz/deta)^2)*D2n

a = Grid.H * Grid.y_i / (Grid.H - 2.0 * Grid.y_i);
b = 1.0 + 2.0 * a / Grid.H;

dzdy   = zeros(Grid.ny, 1);
d2zdy2 = zeros(Grid.ny, 1);
for k = 1:Grid.ny
    dzdy(k)   =  a * (1 + b) / (a + Grid.y(k))^2;
    d2zdy2(k) = -2 * a * (1 + b) / (a + Grid.y(k))^3;
end
Grid.dzdy   = dzdy;
Grid.d2zdy2 = d2zdy2;

Dn  = -chebDiff(Grid.ny);
D2n = Dn * Dn;


% Grid.D1 = diag(dzdy)   * Dn;
% Grid.D2 = diag(d2zdy2) * Dn + diag(dzdy.^2) * D2n;

%Antonis 2026%
Grid.D1 =  Dn;
Grid.D2 =  D2n;



end
