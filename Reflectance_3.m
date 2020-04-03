clear all;

%Constants
h = 6.626e-34; %Plank constant [J.s]
hbar = h/(2*pi);
c = 299792458; % light speed [m/s]
Ta = 300; %Room temperature
%Initial parameters
%cav = 0; % Cavity architecture - 0 for lower refractive index and 1 for a high refractive index
c1 = 0; %x=1 - Al mole fraction in first layer (0 - 1) - in air contact!
c2 = 0.9; %x=0.2 - Al mole fraction in second layer (0 - 1)
c3 = 0; %x=0.7 - Al mole fraction in the cavity (0 - 1)
c4 = 0; % x=0 - Al mole fraction in the cap layer, substrate (0 - 1)
cqw = 0.13; %In mole fraction in InGaAs for QW refractive index; 0 is the GaAs.
n0 = 1.45; %Refractive index of the external medium
ncs = 20; %Number of pair of layers for the upper DBR 
nci = 20; %Number of pair of layers for the bottom DBR 
lambdaR1 = 898; %Ressonant wavelength at 300K [nm]
s = 1; %Cavity order (multiple of lambda/2)
N=3; % QWs number
xqw = 7.5; %QWs thickness [nm] at 300K (all will be equal)
xb = 5; %Barreir thickness [nm] at 300K (all will be equal)
V = 0; %Aplied bias (0 - 10V)
T = 10; %Sample temperature
teta0 = 0*(pi/180); %Light incident angle [rad], in relation to normal of the sample
phi = 0*(pi/180); %Polarization angle (0 = TE, pi/2 = TM) [rad]
shift = 0; %Cavity shift [nm]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change on the ressonance wavelength as function of applied bias
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Due to applied bias, the new ressonance wavelength position at 300K is [nm]:');
lambdaR = lambdaR1 + 4.8825E-2*V - 1.7137E-2*V^2 %Comprimento de onda ressonante - [nm]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Refractive index for the new resonance (due to applied bias)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Refractive index for ressonant wavelength:');
disp('1st DBR layer:');
n1R = Dispersion_2(c1,lambdaR,T) %Dispersao da primeira camada
disp('2nd DBR layer:');
n2R = Dispersion_2(c2,lambdaR,T) %Dispersao da segunda camada
disp('Cavity and barreirs:');
n3R = Dispersion_2(c3,lambdaR,T) %Cavity
disp('Caplayer and Substrate:');
n4R = Dispersion_2(c4,lambdaR,T) %Substrate, QW, Caplayer
disp('QW:');
[nqwR lambInGaAs EgInGaAs EgGaAs] = InGaAsDispersion(cqw,lambdaR,T);
nqwR % In mole fraction in the cap layer, substrate and QWs (0 - 1)

disp('The InGaAs wavelength gap is:');
lambInGaAs
disp('The InGaAs Gap Energy [eV] at defined temperature is:');
EgInGaAs
disp('The GaAs Gap Energy [eV] at defined temperature is:');
EgGaAs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Dispersion curves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda_1 = linspace(700, lambdaR - 30, 500); %Spectrum
lambda_2 = linspace(lambdaR - 30, lambdaR + 30, 5000); %Spectrum
lambda_3 = linspace(lambdaR + 30, 1100, 500); %Spectrum
lambda = [lambda_1 lambda_2 lambda_3];

for b = 1 : length(lambda)
   n_1(b) = Dispersion_2(c1,lambda(b),T); %Dispersao da primeira camada 
   n_2(b) = Dispersion_2(c2,lambda(b),T); %Dispersao da segunda camada
   n_3(b) = Dispersion_2(c3,lambda(b),T); %Dispersao da cavidade 
   n_4(b) = Dispersion_2(c4,lambda(b),T); %Dispersao do substrato
   [n_qw(b) lambInGaAs EgInGaAs EgGaAs] = InGaAsDispersion(cqw,lambda(b),T); %QW Dispersion
end

figure(1)
plot(lambda,n_1,':',lambda,n_2,'--',lambda,n_3,'r',lambda,n_4,'k',lambda,n_qw,'--k' )
xlabel('Wavelength (nm)')
ylabel('Refractive Index')
title('Dispersion for Al_{x}Ga_{x-1}As','FontSize',12)
legend('1st layer','2nd layer','Cavity','GaAs','InGaAs')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%QW absorption coefficient
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Cavity absorption coefficient [m-1] and gap wavelength (nm) at defined temperature:');
[abs lamb_gap4] = Absorption(c3,lambdaR,T) %Absorption coefficien, Gap wavelength

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Building the sample
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xT xqwT xbT xcavHT xcavT Lc] = Sample_2(s,ncs,nci,N,xqw,xb,c1,c2,c3,c4,cqw,lambdaR,T,teta0,shift);
 
ac=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Manual adjust of the layers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x1RT = 64.7;
% x1 = round(Expansion(s,x1RT,c4,Ta,T),1) %New thickness of caplayer
% x2RT = 78.2;
% x2 = round(Expansion(s,x2RT,c4,Ta,T),1) %New thickness of caplayer
% ac = round(-78.2*0.05,1);
% xcavHT = 0;
% xcavRT = 55.5;
% xcavT = round(Expansion(s,xcavRT,c4,Ta,T),1); %New thickness of caplayer
% xqwRT = 7.5;
% xqwT = round(Expansion(s,xqwRT,c4,Ta,T),1); %New thickness of caplayer
% xbRT = 5;
% xbT = round(Expansion(s,xbRT,c4,Ta,T),1); %New thickness of caplayer
% xcl = 0;
% xclT =xcl;
% 
% xT = [x1 x2+ac x1 x2+ac x1 x2+ac x1 x2+ac x1 x2+ac x1 x2+ac x1 x2+ac x1 x2+ac x1 x2+ac x1 x2+ac x1 x2+ac x1 x2+ac x1 x2+ac x1 x2+ac x1 x2+ac x1 x2+ac xbT xqwT xbT xqwT xbT xqwT xcavT xqwT xbT xqwT xbT xqwT xbT x2 x1 x2 x1 x2 x1 x2 x1 x2 x1 x2 x1 x2 x1 x2 x1 x2 x1 x2 x1 x2 x1 x2 x1 x2 x1 x2 x1 x2 x1 x2 x1 x2 x1 x2 x1 x2]; %Manual adjust of the DBR thickness
% 
% Lc = 2*(3*xqwT + 3*xbT) + xcavT;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ncs == 0 %DBR only
    if n1R == n4R 
        xcl = 0;
        x1 = xT(1); % DBR 1st layer
        x2 = xT(2); % DBR 2nd layer
        disp('Do you have a DBR only without caplayer');
    else
        xcl = xT(1); % Caplayer
        x1 = xT(2); % DBR 1st layer
        x2 = xT(3); % DBR 2nd layer
        disp('Do you have a cavity with caplayer');
    end
else
    if nci == 0 %DBR only
        if n1R == n4R 
            xcl = 0;
            x1 = xT(1); % DBR 1st layer
            x2 = xT(2); % DBR 2nd layer
            disp('Do you have a DBR only without caplayer');
        else
            xcl = xT(1); % Caplayer
            x1 = xT(2); % DBR 1st layer
            x2 = xT(3); % DBR 2nd layer
            disp('Do you have a cavity with caplayer');
        end
    else %Microcavity
        if xT(1) <= 1
            xcl = xT(1); % Caplayer
            x1 = xT(2); % DBR 1st layer
            x2 = xT(3); % DBR 2nd layer
            disp('Do you have a cavity with caplayer');
        else
            xcl = 0; % Caplayer
            x1 = xT(1); % DBR 1st layer
            x2 = xT(2); % DBR 2nd layer
            disp('Do you have a cavity without caplayer');
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Reflectance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for b=1:length(lambda) %Scanning the spectra
    m = [1 0; 0 1];
    na= n0;
    teta = teta0;

    for j = 1:length(xT)%(2*ncs)+1 %Scanning the structure
        if xT(j) == xcl %Cap layer
            n = n_4(b);
            [abs lamb_gap4] = Absorption(c4,lambda(b),T); %Absorption coefficient [m-1]
            tetacl =  asin((na/n)*sin(teta));
            p = n*cos(tetacl)*cos(phi) + (n/cos(tetacl))*sin(phi);
            kL = ((2*pi*n/lambda(b)) + i*(abs*(1e-9)/2))*xT(j)/cos(tetacl);
            M = [cos(kL*cos(tetacl)) -(i/p)*sin(kL*cos(tetacl)); -i*p*sin(kL*cos(tetacl)) cos(kL*cos(tetacl))];
            m = m*M;
            na = n_4(b);
            teta = tetacl;
        else %Structure
            if xT(j) == x1 %1st DBR layer
                n = n_1(b);
                [abs lamb_gap1] = Absorption(c1,lambda(b),T); %[m-1]
                teta_1 =  asin((na/n)*sin(teta));
                p = n*cos(teta_1)*cos(phi) + (n/cos(teta_1))*sin(phi);
                kL = ((2*pi*n/lambda(b)) + i*(abs*(1e-9)/2))*xT(j)/cos(teta_1);
                M = [cos(kL*cos(teta_1)) -(i/p)*sin(kL*cos(teta_1)); -i*p*sin(kL*cos(teta_1)) cos(kL*cos(teta_1))];
                m = m*M;
                na = n_1(b);
                teta = teta_1;
            end
            if xT(j) == x2 %2nd DBR layer
                n = n_2(b);
                [abs lamb_gap2] = Absorption(c2,lambda(b),T); %[m-1]
                teta_2 =  asin((na/n)*sin(teta));
                p = n*cos(teta_2)*cos(phi) + (n/cos(teta_2))*sin(phi);
                kL = ((2*pi*n/lambda(b)) + i*(abs*(1e-9)/2))*xT(j)/cos(teta_2);
                M = [cos(kL*cos(teta_2)) -(i/p)*sin(kL*cos(teta_2)); -i*p*sin(kL*cos(teta_2)) cos(kL*cos(teta_2))];
                m = m*M;
                na = n_2(b);
                teta = teta_2;
            end
            if ac == 0
                
            else
                if xT(j) == x1+ac %1st DBR layer
                    n = n_1(b);
                    [abs lamb_gap1] = Absorption(c1,lambda(b),T); %[m-1]
                    teta_1 =  asin((na/n)*sin(teta));
                    p = n*cos(teta_1)*cos(phi) + (n/cos(teta_1))*sin(phi);
                    kL = ((2*pi*n/lambda(b)) + i*(abs*(1e-9)/2))*xT(j)/cos(teta_1);
                    M = [cos(kL*cos(teta_1)) -(i/p)*sin(kL*cos(teta_1)); -i*p*sin(kL*cos(teta_1)) cos(kL*cos(teta_1))];
                    m = m*M;
                    na = n_1(b);
                    teta = teta_1;
                end 
                if xT(j) == x2+ac %2nd DBR layer
                    n = n_2(b);
                    [abs lamb_gap2] = Absorption(c2,lambda(b),T); %[m-1]
                    teta_2 =  asin((na/n)*sin(teta));
                    p = n*cos(teta_2)*cos(phi) + (n/cos(teta_2))*sin(phi);
                    kL = ((2*pi*n/lambda(b)) + i*(abs*(1e-9)/2))*xT(j)/cos(teta_2);
                    M = [cos(kL*cos(teta_2)) -(i/p)*sin(kL*cos(teta_2)); -i*p*sin(kL*cos(teta_2)) cos(kL*cos(teta_2))];
                    m = m*M;
                    na = n_2(b);
                    teta = teta_2;
                end
            end
            
            if xT(j) == xcavHT %Cavity layer
                n = n_3(b); 
                [abs lamb_gap3] = Absorption(c3,lambda(b),T); %[m-1]
                teta_cavH =  asin((na/n)*sin(teta));
                p = n*cos(teta_cavH)*cos(phi) + (n/cos(teta_cavH))*sin(phi);
                kL = ((2*pi*n/lambda(b)) + i*(abs*(1e-9)/2))*xT(j)/cos(teta_cavH);
                M = [cos(kL*cos(teta_cavH)) -(i/p)*sin(kL*cos(teta_cavH)); -i*p*sin(kL*cos(teta_cavH)) cos(kL*cos(teta_cavH))];
                m = m*M;
                na = n_3(b);
                teta = teta_cavH;
            end
            if xT(j) == xqwT %QW
                n = n_qw(b); 
                [abs lamb_gap4] = Absorption(c4,800,T); %[m-1] for GaAs
                abs = 6e4;%[m-1] ~6e6 for InGaAs at 10K
                teta_qw = asin((na/n)*sin(teta));
                p = n*cos(teta_qw)*cos(phi) + (n/cos(teta_qw))*sin(phi);
                kL = ((2*pi*n/lambda(b)) + i*(abs*(1e-9)/2))*xT(j)/cos(teta_qw);
                M = [cos(kL*cos(teta_qw)) -(i/p)*sin(kL*cos(teta_qw)); -i*p*sin(kL*cos(teta_qw)) cos(kL*cos(teta_qw))];
                m = m*M;
                na = n_qw(b);
                teta = teta_qw;
            end
            if xT(j) == xbT %Barreir
                n = n_3(b); 
                [abs lamb_gap3] = Absorption(c3,lambda(b),T); %[m-1]
                teta_b =  asin((na/n)*sin(teta));
                p = n*cos(teta_b)*cos(phi) + (n/cos(teta_b))*sin(phi);
                kL = ((2*pi*n/lambda(b)) + i*(abs*(1e-9)/2))*xT(j)/cos(teta_b);
                M = [cos(kL*cos(teta_b)) -(i/p)*sin(kL*cos(teta_b)); -i*p*sin(kL*cos(teta_b)) cos(kL*cos(teta_b))];
                m = m*M;
                na = n_3(b);
                teta = teta_b;
            end
            if xT(j) == xcavT %Cavity layer
            n = n_3(b); %Cavity layer
            [abs lamb_gap3] = Absorption(c3,lambda(b),T); %[m-1]
            teta_cav =  asin((na/n)*sin(teta));
            p = n*cos(teta_cav)*cos(phi) + (n/cos(teta_cav))*sin(phi);
            kL = ((2*pi*n/lambda(b)) + i*(abs*(1e-9)/2))*xT(j)/cos(teta_cav);
            M = [cos(kL*cos(teta_cav)) -(i/p)*sin(kL*cos(teta_cav)); -i*p*sin(kL*cos(teta_cav)) cos(kL*cos(teta_cav))];
            m = m*M;
            na = n_3(b);
            teta = teta_cav;
            end
        end
    end

 A = m(1,1);
 B = m(1,2); 
 C = m(2,1); 
 D = m(2,2);
 ns = Dispersion_2(c4,lambda(b),T);
 p0 = n0*cos(teta0)*cos(phi) + (n0/cos(teta0))*sin(phi);

 if xT(length(xT)) == x1 || xT(length(xT)-1)%Last 2nd DBR layer
    na = Dispersion_2(c1,lambda(b),T);
    teta_s = asin((na/ns)*sin(teta));
 end
 
 if xT(length(xT)) == x2 || xT(length(xT))%Last 2nd DBR layer
    na = Dispersion_2(c2,lambda(b),T);
    teta_s = asin((na/ns)*sin(teta));
 end
 
 ps = ns*cos(teta_s)*cos(phi) + (ns/cos(teta_s))*sin(phi);
  
 r = (A*p0 + B*p0*ps - C - D*ps)/(A*p0 + B*p0*ps + C + D*ps);
 t = 2*p0/(A*p0 + B*p0*ps + C + D*ps);
 R(b) = r*(r');
 Transm(b) = t*(t');
 %absMic(b) = R(b) - Transm(b); %Absorption
 phase(b) = angle(((n0*B) - C)/((n0*B)+C)); %in [rad]
end
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MQW spectral position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lamb_QW = 888.37;%3QW (InGaAs/GaAs) emission at 10K. % (1e9)*h*c/((1.6e-19)*(1.5187+0.0426+0.0137)) %Wavelength of the QW emission by QWS software [nm]
linewidth = 1.2;%FWMH at 10K.% (1e9)*h*c*((1.6e-19)*(0.076-0.0426+0.0194-0.0137))/((1.6e-19)*(1.5187+0.051+0.0174))^2%QW linewidth emission [nm]

G = 3*normpdf(lambda,lamb_QW,linewidth);% The emission peak of the gain medium


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2)
plot(lambda,R,'r', lambda,G,'b');%, lambda,Transm,'k')
xlabel('Wavelength (nm)')
ylabel('Reflectivity')
title('Reflectivity based in Al_{x}Ga_{x-1}As microcavity','FontSize',12)

data = [lambda' R'];
save('Reflectance.dat','data','-ascii')
% figure(4)
% plot(lambda,absMic)
% xlabel('Wavelength (nm)')
% ylabel('Absorption')
% title('Absorption for the microcavity','FontSize',12)

figure(5)
plot(lambda,phase)
xlabel('Wavelength (nm)')
ylabel('Phase')
title('Phase change spectrum','FontSize',12)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%10nm GaAs QW emission wavelength as a function of Al mole fraction in
%Al_xGa(1-x)As barreir
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wavo = 870;
% Aw = -16.92;
% tw = 0.185;
% Bw = -8.2;
% twi = 0.0223;
% 
% disp('The QW emission is about [nm] (300K):');
% lambda_e = wavo + Aw*(1 - exp(-c3/tw)) + Bw*(1 - exp(-c3/twi)) % [nm]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Penetration length of the optical mode and effective refractive index
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if n1R > n2R
    nH = n1R;
    nL = n2R;
    LH = lambdaR/(4*n1R);
    LHT = Expansion(s,LH,c1,Ta,T);%New thickness
    LL = lambdaR/(4*n2R);
    LLT = Expansion(s,LL,c2,Ta,T);%New thickness
else
    nH = n2R;
    nL = n1R;
    LH = lambdaR/(4*n2R);
    LHT = Expansion(s,LH,c2,Ta,T);%New thickness
    LL = lambdaR/(4*n1R);
    LLT = Expansion(s,LL,c1,Ta,T);%New thickness
end
q = n3R/nL;
p = nL/nH;
a_top = nH/n0;
a_bottom = nH/n4R;

a = a_top;
m = 2*ncs;
disp('The penetration length [nm] for to defined temperature is:');
Lp_top = (lambdaR/(4*n3R))*(q/(1-p))*(1+(a^2)*(p^(m-1))*(1-(p^m)))/(1+((q^2)*(a^2)*(p^(2*m-2))))

a = a_bottom;
m = (2*nci) + 1;
Lp_bottom = (lambdaR/(4*n3R))*(q/(1-p))*(1+(a^2)*(p^(m-1))*(1-(p^m)))/(1+((q^2)*(a^2)*(p^(2*m-2))))

nDBR = ((nH*LHT)+(nL*LLT))/(LHT + LLT);

disp('The effective refractive index to the defined temperature is:');
neff = (nDBR*(Lp_top + Lp_bottom) + lambdaR)/(Lp_top + Lp_bottom + (lambdaR/n3R))

disp('The cavity couplilng coefficient to the defined temperature is:');
Dn = sqrt((n1R - n2R)^2);
kapa= pi*Dn/(2*lambdaR)

disp('The cavity FSR [nm] to the defined temperature is:');
FSR = (lambdaR^2)/(n3R*Lc)

% disp('The threshold gain for the defined temperature is:');
% alfa_a = ;%loss coefficient in the active region
% C = ; %power coupling efficiency between the active and passive sections
% gth = alfa_a+(1/xqwT)*(xcavT*alfa_a - ln(max(R)) +ln(1/C))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Cavity resonance, linewidth and Q factor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambmin = 900-20;
lambmax = 900+20

u=1;
for b=1:length(lambda)
    if lambda(b) >= lambmin
        q = b; %marca o o contador do primeiro termo da varredura
        break
    end
end

for b=1:length(lambda)
    if lambda(b) >= lambmax
        ss = b; %marca o o contador do ultimo termo da varredura
        break
    end
end

for ff = 1 : (ss-q+1) %Tomando o intervalo de varredura 
    lambv(ff) = lambda(q); %Wavelength range of the interval
    Reflec(ff) = R(q); %Reflectance range of the interval
    q = q+1;
end

 figure(6)
 plot(lambv,Reflec)

for b=1:length(lambda)
    if lambda(b) >= lambmin
        q = b; %marca o o contador do primeiro termo da varredura
        break
    end
end

for b=1:length(lambda)
    if lambda(b) >= lambmax
        s = b; %marca o o contador do ultimo termo da varredura
        break
    end
end

for b=q : s %Encontrando o contador da ressonancia
    if R(b) == min(Reflec) 
       u = b;
        break
    end
end

Rg = max (R);
Rp = R(u);
Rm = Rp + ((Rg - Rp)/2);

for b = q : u %Encontrando os contadores da FWHM
    if R(b) >= Rm
        yg = b;
    else
        yp = b;
        break
    end
end

for b = u : s %Encontrando os contadores da FWHM
    if R(b) <= Rm
        wp = b;
    else
        wg = b;
        break
    end
end

%Cavity Figure of Merit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Spontaneus emission coupling factor - beta - Y. Yamamoto, S. Machida, Y. Horikoshi, K. Igeta and G. Bjork, OPTICS COMMUNICATIONS, 80, 5(6), 1991. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

radii = (20e-6)/2; %Pump radii spot size [nm]
meff = (n1R+n2R)/((n1R-n2R)*2);
ne = 2*n1R*n2R/(n1R+n2R);
disp('Effective cavity length [nm] to the defined temperature is:');
Leff = n3R*Lc + 2*ne*meff*(xT(2) + xT(3)) %Cavity effective length [nm]
%Veff = Leff*(pi*radii^2) %[nm^3]
Veff = Lc*(1e-9)*(pi*radii^2) %[m^3]
%Veff = 2.55*(lambdaR/(2*n3R))^3 %[nm^3] - http://optoelectronics.eecs.berkeley.edu/ey1998iee1456.pdf
Dlambda = 10e-9;  %spontaneous emission linewidth [m]
beta = ((lambdaR*(1e-9))^4)/(4*(pi^2)*Veff*Dlambda*n3R)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Cavity linewidth - FHWM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('The cavity FHWM [nm] is:')
FHWM = ((lambda(wg)+ lambda(wp)) - (lambda(yg) + lambda(yp)))/2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Cavity resonance position - Resonance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('The wavelength cavity resonance [nm] is:')
Resonance = lambda(u)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Cavity quality factor - Q
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('The cavity quality factor (Q) is:')
Q = lambda(u)/FHWM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The stop band width
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('The stop band width [nm] for the defined temperature is about:');
StopBand =4*Resonance*sqrt((n1R - n2R)^2)/(pi*(n1R + n2R))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The Purcell factor - https://arxiv.org/pdf/1107.4601.pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('The Purcell Factor for the defined temperature is about:');
Fp = (3/(4*(pi^2)))*Q*((Resonance*(1e-9)/n3R)^3)/Veff
Fp = Q*((Resonance*(1e-9)/n3R)^3)/(4*pi*Veff) %Optical and Quantum Electronics, 24, S245 (1992).
tau_o = 1e-9;
tau = Fp*tau_o; %Spntaneous emission rate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Cavity lifetime
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('The cavity lifetime for the defined temperature is about:');
t_rt=(2*n3R*Lc*(1e-9))/c; %Cavity time round trip
t_cav = t_rt/max(R)
t_cav = ((Resonance*(1e-9))^2)/(2*pi*c*FHWM*(1e-9))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Micropilar energy emission - T. Jakubczyk et. al.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Micropilar energy emission for the defined temperature is about:');
eo = 8.85418782e-12; %Vaccum permitivity
Radii = 5e-6; %Radii micropilar [m]
e = n3R^2; % Cavity Relative Permitivity
Eo = h*(c/(Resonance*1e-9)); %Matching! Emission with the same energy from the cavity [J]
EP = sqrt((Eo^2) + ((hbar*c)^2)/(e*Radii^2)); %Micropilar energy [J]
lambdaP = (1e9)*(c*h/EP) %Micropilar wavelength [nm]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Cavity coupling constant - IEEE J. SELECTED TOPICS IN QUANTUM ELECTRONICS, 21(6), 1503209(2015)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kc = c*(1-max(R))/(2*n3R*Lc*sqrt(max(R)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Threshold inversion density - Lasers, Siegman (pg. 41)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('The threshold inversion photon density for the defined temperature is about:');
tau_sp = 0.5e-9;%spontaneous emission lifetime [s]
g_rad = 1/tau_sp; %radiative decy rate [Hz]
DN = ((2*(pi^2)*c*FHWM*(1e-9)/((Resonance*(1e-9))^2))/(((Resonance*(1e-9)/n3R)^2)*g_rad*Lc))*log(1/max(R))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Total mirror loss - http://www.hunter.cuny.edu/physics/courses/physics852/repository/files/8%20%20Lasers.pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('The cavity loss [nm-1] for the defined temperature is about:');
Cav_loss = (1/(2*Lc))*(log(1/(max(R))^2)) %[nm-1]


Dist = sqrt((lambda(u) - lamb_QW)^2);

    z1 = beta %MAXIMIZE
    z2 = Q %MAXIMIZE
    z3 = Fp %MAXIMIZE
    z4 = DN %MINIMIZE
    z5 = Cav_loss %MINIMIZE
    z6 = Dist %MINIMIZE

