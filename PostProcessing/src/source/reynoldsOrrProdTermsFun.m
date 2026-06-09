
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                      Jordi Casacuberta                               %%
%%                            2021                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Function to compute the production term of the Reynolds-Orr equation
%
% > Inputs: base flow (sBF), perturbation field (foZu/v/w)
% > Outputs: production (sub-)terms (sROProd)
% 

%% 0.- Data loading and definition of inputs 

%% Modal information
%Number of modes
aux.Nh   = 1;
%Primary spanwise wavenumber
aux.beta = 9.64; 

%% Base flow
sBF.u  = BF.u;
sBF.v  = BF.v;
sBF.w  = 0.*BF.u;
sBF.ux = BF.dxU;
sBF.uy = BF.dyU;
sBF.vx = BF.dxV;
sBF.vy = BF.dyV;
sBF.wx = 0.*BF.dxU;
sBF.wy = 0.*BF.dyU;

%% Perturbation amplitude functions
foZu.A(:,:,1)     = abs(u);
foZu.Phase(:,:,1) = angle(u);
foZv.A(:,:,1)     = abs(v);
foZv.Phase(:,:,1) = angle(v);
foZw.A(:,:,1)     = 0.*abs(u);
foZw.Phase(:,:,1) = 0.*angle(u);

%% 1.- Computation of Reynolds-Orr (global) production term (inherently integrated in z)

%Total production term (iterate along Fourier modes)
for k = 1:aux.Nh
    
    %Total production term (integrated in z)
    sROProd.P(:,:,k)  = -  ( 2 * real(foZu.A(:,:,k) ./ 2 .* exp(1i * foZu.Phase(:,:,k)) .* foZu.A(:,:,k) ./ 2 .* exp(-1i * foZu.Phase(:,:,k))) .* sBF.ux ...
                           + 2 * real(foZu.A(:,:,k) ./ 2 .* exp(1i * foZu.Phase(:,:,k)) .* foZv.A(:,:,k) ./ 2 .* exp(-1i * foZv.Phase(:,:,k))) .* sBF.uy ...
                           + 2 * real(foZv.A(:,:,k) ./ 2 .* exp(1i * foZv.Phase(:,:,k)) .* foZu.A(:,:,k) ./ 2 .* exp(-1i * foZu.Phase(:,:,k))) .* sBF.vx ...
                           + 2 * real(foZv.A(:,:,k) ./ 2 .* exp(1i * foZv.Phase(:,:,k)) .* foZv.A(:,:,k) ./ 2 .* exp(-1i * foZv.Phase(:,:,k))) .* sBF.vy ...
                           + 2 * real(foZw.A(:,:,k) ./ 2 .* exp(1i * foZw.Phase(:,:,k)) .* foZu.A(:,:,k) ./ 2 .* exp(-1i * foZu.Phase(:,:,k))) .* sBF.wx ...  
                           + 2 * real(foZw.A(:,:,k) ./ 2 .* exp(1i * foZw.Phase(:,:,k)) .* foZv.A(:,:,k) ./ 2 .* exp(-1i * foZv.Phase(:,:,k))) .* sBF.wy ...     
                           ) ;
    
    %Add the constant factor (coming from integrating int(1 dz) = lambda_z)
    sROProd.P(:,:,k) = sROProd.P(:,:,k) .* (2*pi/aux.beta);
                                            
    %Fill with 0s inside the solid boundary 
    %if rough.include == true && rough.type(1) == 1
    %    sROProd.P(rough.hPos+1:end,rough.xPos:end,k) = 0;   
    %end
end

%% 4.- Computation of Reynolds-Orr production term decomposed following the
%%     tangential-to-normal decomposition 
for k = 1:aux.Nh
    
    %% Define auxiliary variables 
    
    %Note that the amplitude functions are divided by a factor 2 to
    %account for the original -N to N Fourier series
    
    %Module of base-flow vector field
    aux.moduleBF        = sqrt(sBF.u.^2 + sBF.v.^2 + sBF.w.^2);
    
    %Compute auxiliary functions gamma(x,y) and psi(x,y) and xi(x,y)
    aux.gamma1(:,:,k)   = foZu.A(:,:,k) ./ 2 .* sBF.u .* cos(foZu.Phase(:,:,k)) +...
                          foZv.A(:,:,k) ./ 2 .* sBF.v .* cos(foZv.Phase(:,:,k)) +...
                          foZw.A(:,:,k) ./ 2 .* sBF.w .* cos(foZw.Phase(:,:,k));

    aux.gamma2(:,:,k)   = -foZu.A(:,:,k) ./ 2 .* sBF.u .* sin(foZu.Phase(:,:,k)) +...
                          -foZv.A(:,:,k) ./ 2 .* sBF.v .* sin(foZv.Phase(:,:,k)) +...
                          -foZw.A(:,:,k) ./ 2 .* sBF.w .* sin(foZw.Phase(:,:,k));

    aux.psi1(:,:,k)     = (foZu.A(:,:,k) ./ 2).^2 .* cos(2 .* foZu.Phase(:,:,k)) +...
                          (foZv.A(:,:,k) ./ 2).^2 .* cos(2 .* foZv.Phase(:,:,k)) +...
                          (foZw.A(:,:,k) ./ 2).^2 .* cos(2 .* foZw.Phase(:,:,k));

    aux.psi2(:,:,k)     = -(foZu.A(:,:,k) ./ 2).^2 .* sin(2 .* foZu.Phase(:,:,k)) +...
                          -(foZv.A(:,:,k) ./ 2).^2 .* sin(2 .* foZv.Phase(:,:,k)) +...
                          -(foZw.A(:,:,k) ./ 2).^2 .* sin(2 .* foZw.Phase(:,:,k));

    aux.xi(:,:,k)       = sqrt(aux.gamma1(:,:,k).^2 + aux.gamma2(:,:,k).^2) ./ aux.moduleBF.^2;
    
    %% Elements of the tangential vector
    
    %'Shape' of vector components
    %First component
    aux.shapeTangA(:,:,k) = (sqrt(aux.gamma1(:,:,k).^2 + aux.gamma2(:,:,k).^2)) ./ (aux.moduleBF).^2 .* sBF.u;
    
    %Second component
    aux.shapeTangB(:,:,k) = (sqrt(aux.gamma1(:,:,k).^2 + aux.gamma2(:,:,k).^2)) ./ (aux.moduleBF).^2 .* sBF.v;
    
    %Third component
    aux.shapeTangC(:,:,k) = (sqrt(aux.gamma1(:,:,k).^2 + aux.gamma2(:,:,k).^2)) ./ (aux.moduleBF).^2 .* sBF.w;
    
    %Phase of tangential component
    hlp.phaseTang(:,:,k) = atan2(-aux.gamma2(:,:,k),aux.gamma1(:,:,k));
    
    %Fourier coefficient of vector components
    %First component
    aux.FCTangA(:,:,k) = aux.shapeTangA(:,:,k) .* exp(1i * hlp.phaseTang(:,:,k));
    
    %Second component
    aux.FCTangB(:,:,k) = aux.shapeTangB(:,:,k) .* exp(1i * hlp.phaseTang(:,:,k));
    
    %Third component
    aux.FCTangC(:,:,k) = aux.shapeTangC(:,:,k) .* exp(1i * hlp.phaseTang(:,:,k));
    
    %% Elements of the normal vector
    
    %Phase of vector components
    %First component
    aux.phaseNormA(:,:,k) = atan2(foZu.A(:,:,k)./2 .* sin(foZu.Phase(:,:,k)) - sBF.u .* aux.xi(:,:,k) .* sin(hlp.phaseTang(:,:,k)),...
                                  foZu.A(:,:,k)./2 .* cos(foZu.Phase(:,:,k)) - sBF.u .* aux.xi(:,:,k) .* cos(hlp.phaseTang(:,:,k)));
                              
    %Second component
    hlp.phaseNormB(:,:,k) = atan2(foZv.A(:,:,k)./2 .* sin(foZv.Phase(:,:,k)) - sBF.v .* aux.xi(:,:,k) .* sin(hlp.phaseTang(:,:,k)),...
                                  foZv.A(:,:,k)./2 .* cos(foZv.Phase(:,:,k)) - sBF.v .* aux.xi(:,:,k) .* cos(hlp.phaseTang(:,:,k)));
                              
    %Third component
    hlp.phaseNormC(:,:,k) = atan2(foZw.A(:,:,k)./2 .* sin(foZw.Phase(:,:,k)) - sBF.w .* aux.xi(:,:,k) .* sin(hlp.phaseTang(:,:,k)),...
                                  foZw.A(:,:,k)./2 .* cos(foZw.Phase(:,:,k)) - sBF.w .* aux.xi(:,:,k) .* cos(hlp.phaseTang(:,:,k)));
    
    %'Shape' of vector components
    %First component
    aux.shapeNormA(:,:,k) = sqrt((foZu.A(:,:,k)./2).^2 + sBF.u.^2 .* aux.xi(:,:,k).^2 -2.*sBF.u.*aux.xi(:,:,k).*abs(foZu.A(:,:,k)./2).*cos(foZu.Phase(:,:,k) - hlp.phaseTang(:,:,k)));
    
    %Second component
    hlp.shapeNormB(:,:,k) = sqrt((foZv.A(:,:,k)./2).^2 + sBF.v.^2 .* aux.xi(:,:,k).^2 -2.*sBF.v.*aux.xi(:,:,k).*abs(foZv.A(:,:,k)./2).*cos(foZv.Phase(:,:,k) - hlp.phaseTang(:,:,k)));
    
    %Third component
    hlp.shapeNormC(:,:,k) = sqrt((foZw.A(:,:,k)./2).^2 + sBF.w.^2 .* aux.xi(:,:,k).^2 -2.*sBF.w.*aux.xi(:,:,k).*abs(foZw.A(:,:,k)./2).*cos(foZw.Phase(:,:,k) - hlp.phaseTang(:,:,k)));
    
    %Fourier coefficient of vector components
    %First component
    aux.FCNormA(:,:,k) = aux.shapeNormA(:,:,k) .* exp(1i * aux.phaseNormA(:,:,k));
    
    %Second component
    aux.FCNormB(:,:,k) = hlp.shapeNormB(:,:,k) .* exp(1i * hlp.phaseNormB(:,:,k));
    
    %Third component
    aux.FCNormC(:,:,k) = hlp.shapeNormC(:,:,k) .* exp(1i * hlp.phaseNormC(:,:,k));
    
    %% Modulus of the full tangential and normal vector
    
    %Norm of total vector, i.e., ||v'||(x,y)
    aux.shapeTot(:,:,k) = sqrt((foZu.A(:,:,k)./2).^2 + (foZv.A(:,:,k)./2).^2 + (foZw.A(:,:,k)./2).^2);

    %Norm of tangential vector, i.e., ||vt'||(x,y)
    hlp.shapeTang(:,:,k) = (sqrt(aux.gamma1(:,:,k).^2 + aux.gamma2(:,:,k).^2)) ./ (aux.moduleBF);

    %Norm of normal vector, i.e., ||vn'||(x,y). Pitagorian
    %expression is applicable in complex vector space as well
    hlp.shapeNorm(:,:,k) = sqrt(aux.shapeTot(:,:,k).^2 - hlp.shapeTang(:,:,k).^2);
    
    %% Compute decomposed production terms
    
    %Terms of I1
    aux.I1a(:,:,k) = 2 .* sBF.ux .* real(aux.FCNormA(:,:,k) .* conj(aux.FCNormA(:,:,k)));
    aux.I1b(:,:,k) = 2 .* sBF.uy .* real(aux.FCNormA(:,:,k) .* conj(aux.FCNormB(:,:,k)));
    aux.I1c(:,:,k) = 2 .* sBF.vx .* real(aux.FCNormB(:,:,k) .* conj(aux.FCNormA(:,:,k)));
    aux.I1d(:,:,k) = 2 .* sBF.vy .* real(aux.FCNormB(:,:,k) .* conj(aux.FCNormB(:,:,k)));
    aux.I1e(:,:,k) = 2 .* sBF.wx .* real(aux.FCNormC(:,:,k) .* conj(aux.FCNormA(:,:,k)));
    aux.I1f(:,:,k) = 2 .* sBF.wy .* real(aux.FCNormC(:,:,k) .* conj(aux.FCNormB(:,:,k)));
    
    %I1 itself
    aux.I1(:,:,k)     = - (2*pi/aux.beta) .* (aux.I1a(:,:,k) + aux.I1b(:,:,k) + aux.I1c(:,:,k) + aux.I1d(:,:,k) + aux.I1e(:,:,k) + aux.I1f(:,:,k));
    sROProd.I1(:,:,k) = aux.I1(:,:,k);
    
    %Terms of I2 (hlp contains each term with factor for comparison with In)
    aux.I2a(:,:,k) = 2 .* sBF.ux .* real(aux.FCTangA(:,:,k) .* conj(aux.FCNormA(:,:,k))); hlp.I2a(:,:,k) = - (2*pi/aux.beta) .* aux.I2a(:,:,k);
    aux.I2b(:,:,k) = 2 .* sBF.uy .* real(aux.FCTangA(:,:,k) .* conj(aux.FCNormB(:,:,k))); hlp.I2b(:,:,k) = - (2*pi/aux.beta) .* aux.I2b(:,:,k);
    aux.I2c(:,:,k) = 2 .* sBF.vx .* real(aux.FCTangB(:,:,k) .* conj(aux.FCNormA(:,:,k))); hlp.I2c(:,:,k) = - (2*pi/aux.beta) .* aux.I2c(:,:,k);
    aux.I2d(:,:,k) = 2 .* sBF.vy .* real(aux.FCTangB(:,:,k) .* conj(aux.FCNormB(:,:,k))); hlp.I2d(:,:,k) = - (2*pi/aux.beta) .* aux.I2d(:,:,k);
    aux.I2e(:,:,k) = 2 .* sBF.wx .* real(aux.FCTangC(:,:,k) .* conj(aux.FCNormA(:,:,k))); hlp.I2e(:,:,k) = - (2*pi/aux.beta) .* aux.I2e(:,:,k);
    aux.I2f(:,:,k) = 2 .* sBF.wy .* real(aux.FCTangC(:,:,k) .* conj(aux.FCNormB(:,:,k))); hlp.I2f(:,:,k) = - (2*pi/aux.beta) .* aux.I2f(:,:,k);
    
    %I2 itself and sub-fields
    aux.I2(:,:,k)      = - (2*pi/aux.beta) .* (aux.I2a(:,:,k) + aux.I2b(:,:,k) + aux.I2c(:,:,k) + aux.I2d(:,:,k) + aux.I2e(:,:,k) + aux.I2f(:,:,k));
    sROProd.I2(:,:,k)  = aux.I2(:,:,k);
    
    sROProd.I2a(:,:,k) = - (2*pi/aux.beta) .* (aux.I2a(:,:,k));  sROProd.I2b(:,:,k) = - (2*pi/aux.beta) .* (aux.I2b(:,:,k));
    sROProd.I2c(:,:,k) = - (2*pi/aux.beta) .* (aux.I2c(:,:,k));  sROProd.I2d(:,:,k) = - (2*pi/aux.beta) .* (aux.I2d(:,:,k));
    sROProd.I2e(:,:,k) = - (2*pi/aux.beta) .* (aux.I2e(:,:,k));  sROProd.I2f(:,:,k) = - (2*pi/aux.beta) .* (aux.I2f(:,:,k));
    
    %Terms of I3
    aux.I3a(:,:,k) = 2 .* sBF.ux .* real(aux.FCNormA(:,:,k) .* conj(aux.FCTangA(:,:,k)));
    aux.I3b(:,:,k) = 2 .* sBF.uy .* real(aux.FCNormA(:,:,k) .* conj(aux.FCTangB(:,:,k)));
    aux.I3c(:,:,k) = 2 .* sBF.vx .* real(aux.FCNormB(:,:,k) .* conj(aux.FCTangA(:,:,k)));
    aux.I3d(:,:,k) = 2 .* sBF.vy .* real(aux.FCNormB(:,:,k) .* conj(aux.FCTangB(:,:,k)));
    aux.I3e(:,:,k) = 2 .* sBF.wx .* real(aux.FCNormC(:,:,k) .* conj(aux.FCTangA(:,:,k)));
    aux.I3f(:,:,k) = 2 .* sBF.wy .* real(aux.FCNormC(:,:,k) .* conj(aux.FCTangB(:,:,k)));
    
    %I3 itself
    aux.I3(:,:,k)     = - (2*pi/aux.beta) .* (aux.I3a(:,:,k) + aux.I3b(:,:,k) + aux.I3c(:,:,k) + aux.I3d(:,:,k) + aux.I3e(:,:,k) + aux.I3f(:,:,k));
    sROProd.I3(:,:,k) = aux.I3(:,:,k);
    
    %Terms of I4
    aux.I4a(:,:,k) = 2 .* sBF.ux .* real(aux.FCTangA(:,:,k) .* conj(aux.FCTangA(:,:,k)));
    aux.I4b(:,:,k) = 2 .* sBF.uy .* real(aux.FCTangA(:,:,k) .* conj(aux.FCTangB(:,:,k)));
    aux.I4c(:,:,k) = 2 .* sBF.vx .* real(aux.FCTangB(:,:,k) .* conj(aux.FCTangA(:,:,k)));
    aux.I4d(:,:,k) = 2 .* sBF.vy .* real(aux.FCTangB(:,:,k) .* conj(aux.FCTangB(:,:,k)));
    aux.I4e(:,:,k) = 2 .* sBF.wx .* real(aux.FCTangC(:,:,k) .* conj(aux.FCTangA(:,:,k)));
    aux.I4f(:,:,k) = 2 .* sBF.wy .* real(aux.FCTangC(:,:,k) .* conj(aux.FCTangB(:,:,k)));
    
    %I4 itself and sub-fields
    aux.I4(:,:,k)     = - (2*pi/aux.beta) .* (aux.I4a(:,:,k) + aux.I4b(:,:,k) + aux.I4c(:,:,k) + aux.I4d(:,:,k) + aux.I4e(:,:,k) + aux.I4f(:,:,k));
    sROProd.I4(:,:,k) = aux.I4(:,:,k);
    
    sROProd.I4a(:,:,k) = - (2*pi/aux.beta) .* (aux.I4a(:,:,k));  sROProd.I4b(:,:,k) = - (2*pi/aux.beta) .* (aux.I4b(:,:,k));
    sROProd.I4c(:,:,k) = - (2*pi/aux.beta) .* (aux.I4c(:,:,k));  sROProd.I4d(:,:,k) = - (2*pi/aux.beta) .* (aux.I4d(:,:,k));
    sROProd.I4e(:,:,k) = - (2*pi/aux.beta) .* (aux.I4e(:,:,k));  sROProd.I4f(:,:,k) = - (2*pi/aux.beta) .* (aux.I4f(:,:,k));
    
    %Fill with 0s inside the solid boundary 
    %if rough.include == true && rough.type(1) == 1
    %    sROProd.I1(rough.hPos+1:end,rough.xPos:end,k) = 0;   
    %    sROProd.I2(rough.hPos+1:end,rough.xPos:end,k) = 0; 
    %    sROProd.I3(rough.hPos+1:end,rough.xPos:end,k) = 0; 
    %    sROProd.I4(rough.hPos+1:end,rough.xPos:end,k) = 0; 
    %end
                                            
end

%% Perform volume integral of Reynolds-Orr terms (along x and y)

%Define integration limits
lim.xMin     = 177.6149; %177.6149;
lim.xMax     = 197.6149; %187.6149;
lim.yHeight  = 6;

if rough.include == true
    lim.yMax = (rough.h + rough.deltaH) + lim.yHeight;
else
    lim.yMax = lim.yHeight;
end

lim.xMinPos = find(squeeze(sBF.X(1,1,:)) > lim.xMin,1,'first');
lim.xMaxPos = find(squeeze(sBF.X(1,1,:)) > lim.xMax,1,'first');
lim.yMaxPos = find(squeeze(sBF.Y(:,1,1)) < lim.yMax,1,'first');

%Perfom volumetric integration (in x and y since the term is arranged
%as already integrated in z)

%Integrate production
for k = 1:aux.Nh
    aux.intVar = sROProd.P(:,:,k);
    if rough.include == true
        sROProd.intP(k)    = -trapz(sBF.Y(lim.yMaxPos:rough.hPos,1,1),trapz(squeeze(sBF.X(1,1,lim.xMinPos:lim.xMaxPos)),...
                              aux.intVar(lim.yMaxPos:rough.hPos,lim.xMinPos:lim.xMaxPos),2));
        sROProd.intAbsP(k) = -trapz(sBF.Y(lim.yMaxPos:rough.hPos,1,1),trapz(squeeze(sBF.X(1,1,lim.xMinPos:lim.xMaxPos)),...
                              abs(aux.intVar(lim.yMaxPos:rough.hPos,lim.xMinPos:lim.xMaxPos)),2));
    else
        sROProd.intP(k)    = -trapz(sBF.Y(lim.yMaxPos:end,1,1),trapz(squeeze(sBF.X(1,1,lim.xMinPos:lim.xMaxPos)),...
                              aux.intVar(lim.yMaxPos:end,lim.xMinPos:lim.xMaxPos),2));
        sROProd.intAbsP(k) = -trapz(sBF.Y(lim.yMaxPos:end,1,1),trapz(squeeze(sBF.X(1,1,lim.xMinPos:lim.xMaxPos)),...
                              abs(aux.intVar(lim.yMaxPos:end,lim.xMinPos:lim.xMaxPos)),2));
    end
end

%Integrate perturbation kinetic energy
for k = 1:aux.Nh
    if rough.include == true
        aux.intTot(k)    = -trapz(sBF.Y(lim.yMaxPos:rough.hPos,1,1),trapz(squeeze(sBF.X(1,1,lim.xMinPos:lim.xMaxPos)),...
                              0.5*aux.shapeTot(lim.yMaxPos:rough.hPos,lim.xMinPos:lim.xMaxPos,k).^2,2)); 
        aux.intTang(k)   = -trapz(sBF.Y(lim.yMaxPos:rough.hPos,1,1),trapz(squeeze(sBF.X(1,1,lim.xMinPos:lim.xMaxPos)),...
                              0.5*hlp.shapeTang(lim.yMaxPos:rough.hPos,lim.xMinPos:lim.xMaxPos,k).^2,2));
        aux.intNorm(k)   = -trapz(sBF.Y(lim.yMaxPos:rough.hPos,1,1),trapz(squeeze(sBF.X(1,1,lim.xMinPos:lim.xMaxPos)),...
                              0.5*hlp.shapeNorm(lim.yMaxPos:rough.hPos,lim.xMinPos:lim.xMaxPos,k).^2,2));
    else
        aux.intTot(k)    = -trapz(sBF.Y(lim.yMaxPos:end,1,1),trapz(squeeze(sBF.X(1,1,lim.xMinPos:lim.xMaxPos)),...
                              0.5*aux.shapeTot(lim.yMaxPos:end,lim.xMinPos:lim.xMaxPos,k).^2,2));
        aux.intTang(k)   = -trapz(sBF.Y(lim.yMaxPos:end,1,1),trapz(squeeze(sBF.X(1,1,lim.xMinPos:lim.xMaxPos)),...
                              0.5*hlp.shapeTang(lim.yMaxPos:end,lim.xMinPos:lim.xMaxPos,k).^2,2));
        aux.intNorm(k)   = -trapz(sBF.Y(lim.yMaxPos:end,1,1),trapz(squeeze(sBF.X(1,1,lim.xMinPos:lim.xMaxPos)),...
                              0.5*hlp.shapeNorm(lim.yMaxPos:end,lim.xMinPos:lim.xMaxPos,k).^2,2));
    end
end

%% Plot results

% %Limits to plot 
% minX       = 175.62;
% maxX       = 190.62;
% maxY       = 4; 
% 
% pathPng = ['InputsOutputs/pngFiles/' generation '/' refinement '/' test '/CFI'];
% 
% if mainFig == true
%     reynoldsOrrProductionPlotterFun(sBF,sROProd,hlp,rough,pathPng,minX,maxX,maxY,refinement);
% end
% 
% clear aux hlp
