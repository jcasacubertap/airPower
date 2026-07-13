function [eigenvalue,index] = EVfilter_channel(EV,V,y,omega,beta,delta99,Mach)
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
eigen   = EV;
for  ii =1:length(eigen)
    if  real(eigen) < 0 
        eigen(ii) = 0 ;
    end     
    if (omega/real(eigen(ii)))>1   % faster waves than flow
        eigen(ii) = 0 ;
    end   
    if real(eigen(ii))==inf   % no infty
        eigen(ii) = 0;
    end
    if (abs(imag(eigen(ii))))>0.3   % unnatural growth rate
        eigen(ii) = 0;
    end
end 

[~,index]=max(real(eigen));   % choose the shortest mode
eigenvalue = eigen(index);

end