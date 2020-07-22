function [nInGaAs lambInGaAs EgInGaAs EgGaAs] = InGaAsDispersion(cIn,lambda,T) %InGaAs QW Refractive index

h = 6.626e-34; %Plank constant [J.s]
v = 299792458; %vacuum light speed [m/s]
e = 1.6e-19; %Elementary electron charge [C]
AInGaAs = 8.95;
BInGaAs = 2.054;
CInGaAs = 624.5;

[a_GaAs lambGap_GaAs] = Absorption(0,lambda,T);

EgGaAs = h*v/(e*lambGap_GaAs*1e-9); %Gap energy for GaAs [eV]

EgInGaAs = EgGaAs - (cIn*1.56) + (cIn^2)*0.494; %InGaAs Gap Energy [eV]

lambInGaAs = h*v*(1e9)/(e*EgInGaAs); %Wavelength InGaAs Gap [nm]

nInGaAs = sqrt(AInGaAs + BInGaAs/(1-(CInGaAs*EgGaAs/(lambda*EgInGaAs))^2)); %3.5484;
end