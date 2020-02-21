%function [Eg Ex] = QWEnergy(x,T) % AlGaAs Gap Energy as function of temperature - Al concentration, Temperature (K)

Ta=300;
x=1;
T=300;
%Constants
e = 1.6e-19;
h = 6.626e-34;
c = 299792458;

%At 300K
if x<=0.45 %Direct gap
        Egoq = (1.6e-19)*(1.519+1.23*x); %[J] - Gap Energy to T=0K
        aq = (1.6e-19)*(1e-4)*(4.9 + 0.7*x + 3.7*(x^2));
        tetaq = 250 + 5*x + 260*(x^2);
        pq = 2.33+1.7*x;
        Eg = ((Egoq - (aq*tetaq*(nthroot((1+((2*Ta/tetaq)^pq)), pq)-1)/2)))/e %Gap Energy [eV]
        Ex = (7.52e-22 + 1.0912e-21*x + 8.768e-22*(x)^2)/e; %1S Exciton binding energy [eV]
else % 0.45 < x < 1 -> Gap Indireto
        Egoq = (1.6e-19)*(1.519*(1-x) + 2.239*x - x*(1-x)*(-0.127+ 1.31*x));
        aq = e*((4.9e-4)*(1-x) + (3.82e-4)*x);
        tetaq = (226)*(1-x) + (250)*x;
        pq =  (2.33)*(1-x) + (2.32)*x;
        Eg = ((Egoq - (aq*tetaq*(nthroot((1+((2*Ta/tetaq)^pq)), pq)-1)/2)))/e %Energia de Gap [eV]
        Ex = (7.52e-22 + 1.0912e-21*x + 8.768e-22*(x)^2)/e; %Energia Ligacao Exciton 1S [eV]
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%10A GaAs QW emission wavelength as a function of Al mole fraction in
%Al_xGa(1-x)As barreir at 300K
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wavo = h*c*(1e9)/(e*Eg);%870;
Aw = -16.92;
tw = 0.185;
Bw = -8.2;
twi = 0.0223;

%disp('The QW emission is about [nm] (300K):');
%lambda_e = wavo + Aw*(1 - exp(-x/tw)) + Bw*(1 - exp(-x/twi)); % [nm]
%E_e = h*c/(e*lambda_e*(1e-9)) %[eV]

%Energy Difference at 300K
DE = (Eg-Ex) - E_e

%For any Temperature
if x<=0.45 %Gap Direto
        Egoq = (1.6e-19)*(1.519+1.23*x); %[J] - Gap Energy to T=0K
        aq = (1.6e-19)*(1e-4)*(4.9 + 0.7*x + 3.7*(x^2));
        tetaq = 250 + 5*x + 260*(x^2);
        pq = 2.33+1.7*x;
        Em = (((Egoq - (aq*tetaq*(nthroot((1+((2*T/tetaq)^pq)), pq)-1)/2)))/e) + DE %Energia de Gap [eV]
else % 0.45 < x < 1 -> Gap Indireto
        Egoq = (1.6e-19)*(1.519*(1-x) + 2.239*x - x*(1-x)*(-0.127+ 1.31*x));
        aq = e*((4.9e-4)*(1-x) + (3.82e-4)*x);
        tetaq = (226)*(1-x) + (250)*x;
        pq =  (2.33)*(1-x) + (2.32)*x;
        Em = (((Egoq - (aq*tetaq*(nthroot((1+((2*T/tetaq)^pq)), pq)-1)/2)))/e) + DE %Energia de Gap [eV]
end