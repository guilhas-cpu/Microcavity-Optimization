function n = Dispersion_2(c,lambda,T)

%Sellmeier coefficients

%%%GaAs
    AGaAs = 3.44;
    BGaAs = 7.54;
    CGaAs = 126424;
%%%AlAs
    AAlAs = 2;
    BAlAs =5.4;
    CAlAs =83441;
    
%Sellmeier Equation
nGaAs = sqrt(AGaAs + ((BGaAs*(lambda^2))/((lambda^2) - CGaAs)));
nAlAs = sqrt(AAlAs + ((BAlAs*(lambda^2))/((lambda^2) - CAlAs)));

%%%%%%%%%%%%%%%%%%%%
%Temperature correction
%%%%%%%%%%%%%%%%%%%%
beta = 1.43e-5; %Thermo-optic coefficient for GaAs by fitting

if c == 0 %GaAs only
    no = nGaAs/exp(beta*300); %Refractive index at T=0K
    nGaAs = no*exp(beta*T); %Refractive index at defined temperature
end

n = c*nAlAs +(1-c)*nGaAs;
end