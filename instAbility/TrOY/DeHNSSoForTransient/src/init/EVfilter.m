function [eigenvalue,index] = EVfilter(EV,V,y,omega,beta,delta99,Mach)
% Function which filters out incorrect eigenvalues from an eigenvalue
% spectrum by checking for an exponential fit
%
% The physical disturbance will decay exponentially along "y". Therefore, 
% each eigenfunction is compared to:
%       phi = exp( - i sqrt(alpha_r^2 + beta^2) y )
% where alpha is the eigenvalue corresponding to each eigenfunction 
% (real component only) and beta the spanwise wavenumber used as input 
% for the OS solver.
%
% It is required that the OS solver used a mapping that includes the full 
% exponential decay of the eigenfunction! Therefore, make sure that y_max 
% of the mapping is sufficiently high.
%
% M. Kotsonis, K. Groot & J.Y. Boersma - 2017

if nargin<7
    Mach = 0;
end

power1 = 6;     % Tweaking factor for the slope error
power2 = 2;     % Tweaking factor for the switch error
treshold = 14;  % logarithm of accuracy of the numerical eigenfunction

use = ~isinf(EV);
EV = EV(use);
V = V(:,use);

NC = length(y);     % length of Chebychev matrix
NE = length(EV);    % length of eigenspectrum

% Parse eigenfunctions V
[Vmax,~] = max(abs(V)); % Store maxima for later use
V = V./repmat(Vmax,NC,1); % Normalize
V(isnan(V))=0+0*sqrt(-1);% replace NaN with zeros

% Combined wavenumber from α and β (Balakumar & Malik 1992, eqs. 3.7, 2.19;
% assumes Pr = 3/4 and Re >> 1).
wavenumber = real(sqrt(EV.^2 + beta^2 - (EV-omega).^2 * Mach^2))';

% Loop over eigenvalues (and functions)
Error = inf(1,NE);
SlopeError = zeros(1,NE);
Switches = zeros(1,NE);
y_top = zeros(1,NE);
y_bulge = zeros(1,NE);
for n = 1:NE
    % Skip clearly unphysical eigenvalues early.
    % The α_r range allows negative values: stationary CFI on a swept BL
    % gives α_r = -β·W_edge/U_edge < 0 when W_edge and β have the same sign.
    % The exponential-decay tests below use wavenumber = sqrt(α_r²+β²) which
    % is sign-independent, so the rest of the filter handles negative α_r
    % correctly.
    if abs(real(EV(n))) > 1.0 || wavenumber(n) < 1e-10
        continue
    end

    %% Limits in y

    % Determine at which y the numerical accuracy of the eigenfunction will no
    % longer produce an exponential trendline
    y_top(n) = treshold/wavenumber(n);

    % Determine the point where the exponential decay starts and the bulge ends
    % (this point is not strictly defined, just an arbitrary height)
    y_bulge(n) = delta99 + 0.15*delta99/wavenumber(n);
    
    % Check the values
    if y_bulge(n) > y(1)/2
        y_bulge(n) = y(1)/2;
    end
    if y_top(n) < y_bulge(n)*2
        y_top(n) = inf;
    end


    %% Fit of the exponential trend
    % Exponential decay should be equal to the wavenumber

    % Selection "s1" selects the freestream
    s1 = max([find(y<=y_top(n),  1),        find(y<=max(y)*0.85,1)]) ...    % Always skip upper 15%
       : min([find(y>=y_bulge(n),1,'last'), NC-1]);                         % Always skip last node
    EF = log(abs(V(s1,n)));
   
    % Per-segment comparison
    SlopeDifference = wavenumber(n) + diff(EF)./diff(y(s1));
    SegmentLength = -diff(y(s1)); % Weighting factor
    SlopeError_Segment = mean(abs(SegmentLength.*SlopeDifference.^2));

    % Linear least-square fit to determine the slope
    EF_fit = polyfit(y(s1),EF,1);
    SlopeError_LS = abs(wavenumber(n)+EF_fit(1));
    
    % Combine segment comparison and least-square fit
    SlopeError(n) = sqrt(SlopeError_Segment.*SlopeError_LS).^(power1);

    %% Number of switches in gradient

    % Selection "s2" selects the bulge
    s2 = find(y>=y_bulge(n),1,'last'):NC;

    % Count switches, reduced by expected values for each component
    S_imag = countSwitches(imag(V(s2,n)),2);
    S_real = countSwitches(real(V(s2,n)),1);
    S_abs  = countSwitches( abs(V(s2,n)),1);

    % The sum of the number of switches should be zero
    Switches(n) = S_imag+S_real+S_abs;

    % Scale the slope error with the number of switches, using an exponent 
    % because of the logarithmic scaling.
    % Times "power" amplifies the error of multiple switches.
    Error(n) = SlopeError(n).*10.^(Switches(n)*power2);
    

end

%% Select eigenvalue

% Minimum error is the eigenvalue to keep
[~,keep] = min(Error);
eigenvalue = EV(keep);

% Select correct index of original EV list
remain = find(use);
index = remain(keep);


end

function switches = countSwitches(data,expected)
% Count the number of 'switches' in the gradient of 'data', reduced by the
% number of 'expected' switches, up to a minimum of zero.

% Sign of difference of the data segments (hence, removes one element)
sn = sign(diff(data));
% Product of succeeding data segments (also removes one element)
count = sn(1:end-1,:).*sn(2:end,:); % -1 if switched, 1 if not
% Count number of switches
switches = sum(count+1==0);
% Reduce number of switches
switches = abs(switches-expected);
% % Apply minimum of zero
% switches(switches<0) = 0;
end
