function [eigenvalue,index] = EVfilter(EV,V,y,omega,beta,delta99,Mach,mode_select)
% Filter an eigenvalue spectrum for the physical mode by an exponential-decay
% fit: the disturbance decays as phi = exp(-i sqrt(alpha_r^2 + beta^2) y), so
% each eigenfunction is compared against that slope. Requires the OS mapping to
% resolve the full decay (y_max sufficiently high).
%
% mode_select (optional, Ma>1 where Mack F and S modes coexist):
%   'auto' (default) — minimum-error eigenvalue
%   'slow'           — slow mode S (phase speed c_r < 1)
%   'fast'           — fast mode F (phase speed c_r > 1)
%
% M. Kotsonis, K. Groot & J.Y. Boersma - 2017

if nargin<7
    Mach = 0;
end
if nargin<8 || isempty(mode_select)
    mode_select = 'auto';
end

power1 = 6;     % slope-error exponent
power2 = 2;     % switch-error exponent
treshold = 14;  % logarithm of eigenfunction accuracy

use = ~isinf(EV);
EV = EV(use);
V = V(:,use);

NC = length(y);
NE = length(EV);

[Vmax,~] = max(abs(V));
V = V./repmat(Vmax,NC,1);
V(isnan(V))=0+0*sqrt(-1);

% Combined wavenumber of alpha and beta, Balakumar & Malik (1992) eqs (3.7),
% (2.19); assumes Pr = 3/4 and Re >> 1.
wavenumber = real(sqrt(EV.^2 + beta^2 - (EV-omega).^2 * Mach^2))';

Error = inf(1,NE);
SlopeError = zeros(1,NE);
Switches = zeros(1,NE);
y_top = zeros(1,NE);
y_bulge = zeros(1,NE);

% Continuous-spectrum branch points (Balakumar & Malik 1992, eq 3.8): reject
% discrete candidates sitting within branch_tol of a branch cut.
branch_tol  = 0.02;
branch_list = omega;                % vorticity/entropy
if Mach > 0
    branch_list = [branch_list, Mach*omega/(Mach+1), Mach*omega/(Mach-1)];
end

% Propagating modes (TS/Mack) sit on the positive branch (alpha_r ~ omega), but
% a (near-)stationary crossflow mode can have alpha_r < 0 when the crossflow
% points the other way; keep both signs there.
stationary = abs(omega) < 1e-6 * max(1, abs(beta));

for n = 1:NE
    if (~stationary && real(EV(n)) < 0) || abs(real(EV(n))) > 1.0 || wavenumber(n) < 1e-10
        continue
    end
    if any(abs(real(EV(n)) - branch_list) < branch_tol*abs(omega)) && ...
       abs(imag(EV(n))) < branch_tol*abs(omega)
        continue
    end

    %% Limits in y
    y_top(n) = treshold/wavenumber(n);
    y_bulge(n) = delta99 + 0.15*delta99/wavenumber(n);
    if y_bulge(n) > y(1)/2
        y_bulge(n) = y(1)/2;
    end
    if y_top(n) < y_bulge(n)*2
        y_top(n) = inf;
    end

    %% Fit of the exponential trend
    % s1 selects the freestream (skip upper 15% and the last node)
    s1 = max([find(y<=y_top(n),  1),        find(y<=max(y)*0.85,1)]) ...
       : min([find(y>=y_bulge(n),1,'last'), NC-1]);
    EF = log(abs(V(s1,n)));

    SlopeDifference = wavenumber(n) + diff(EF)./diff(y(s1));
    SegmentLength = -diff(y(s1));
    SlopeError_Segment = mean(abs(SegmentLength.*SlopeDifference.^2));

    EF_fit = polyfit(y(s1),EF,1);
    SlopeError_LS = abs(wavenumber(n)+EF_fit(1));

    SlopeError(n) = sqrt(SlopeError_Segment.*SlopeError_LS).^(power1);

    %% Number of switches in gradient
    % s2 selects the bulge; the switch count should be zero
    s2 = find(y>=y_bulge(n),1,'last'):NC;
    S_imag = countSwitches(imag(V(s2,n)),2);
    S_real = countSwitches(real(V(s2,n)),1);
    S_abs  = countSwitches( abs(V(s2,n)),1);
    Switches(n) = S_imag+S_real+S_abs;

    Error(n) = SlopeError(n).*10.^(Switches(n)*power2);
end

%% Select eigenvalue

% Optional slow/fast branch restriction (Ma>1). Admissible windows:
%   slow mode (S): omega < Re(alpha) < omega*M/(M-1)   [c_r in (1-1/M, 1)]
%   fast mode (F): omega*M/(M+1) < Re(alpha) < omega   [c_r in (1, 1+1/M)]
switch lower(mode_select)
    case 'slow'
        if Mach > 1
            slow_acf_ub = omega * Mach / (Mach - 1);
            wrong_branch = real(EV) <= omega | real(EV) >= slow_acf_ub;
        else
            wrong_branch = real(EV) <= omega;
        end
        Error(wrong_branch) = inf;
    case 'fast'
        if Mach > 1
            fast_acf_lb = omega * Mach / (Mach + 1);
            wrong_branch = real(EV) >= omega | real(EV) <= fast_acf_lb;
        else
            wrong_branch = real(EV) >= omega;
        end
        Error(wrong_branch) = inf;
    case 'auto'
    otherwise
        warning('EVfilter:UnknownModeSelect', ...
                'mode_select=%s not recognised; using ''auto''', mode_select);
end

[~,keep] = min(Error);
if ~isfinite(Error(keep))
    warning('EVfilter:NoCandidate', ...
            'mode_select=%s left no admissible candidate; falling back to all eigenvalues', mode_select);
    [~,keep] = min(abs(real(EV) - omega));
end
eigenvalue = EV(keep);

remain = find(use);
index = remain(keep);

end

function switches = countSwitches(data,expected)
% Count sign switches in the gradient of 'data', reduced by 'expected'.
sn = sign(diff(data));
count = sn(1:end-1,:).*sn(2:end,:);   % -1 if switched, 1 if not
switches = sum(count+1==0);
switches = abs(switches-expected);
end
