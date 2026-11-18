function w = jacobian_apply(v_real, A, B, N)
%% License (GNU GENERAL PUBLIC LICENSE v3)
%                  Delft Harmonic Navier-Stokes Solver
%     Copyright (C) 2026 P.C. Boldini, Y. Chen & M. Kotsonis
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
% Newton-Jacobian matrix-vector product in real arithmetic.
%
%   w = jacobian_apply(v_real, A, B, N)
%
%   A = ML - Jh and B = -Ja are the paper's J_h and J_a, built from the nlt_jacobian
%   blocks. With vc = vr + 1i*vi this returns
%       [Re(A*vc + B*conj(vc)); Im(A*vc + B*conj(vc))],
%   which is the real Jacobian [Re(A+B), -Im(A-B); Im(A+B), Re(A-B)] applied to
%   [vr; vi]. v_real stacks [real(vc); imag(vc)]; N is the length of vc.

%% Inputs
% Name   size         unit   explanation
% v_real (2N,1)       [-]    stacked real and imaginary parts [real(vc); imag(vc)] of the complex vector vc
% A      (N,N) sparse [-]    J_h block of the Jacobian, ML - Jh
% B      (N,N) sparse [-]    J_a block of the Jacobian, -Ja, acting on the conjugate
% N      (1,1)        [-]    length of the complex vector vc

%% Outputs
% Name size   unit   explanation
% w    (2N,1) [-]    real Jacobian applied to v_real, [Re(A*vc + B*conj(vc)); Im(...)]

vc = v_real(1:N) + 1i * v_real(N+1:2*N);
wc = A * vc + B * conj(vc);
w  = [real(wc); imag(wc)];
end
