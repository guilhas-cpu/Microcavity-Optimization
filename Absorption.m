function [a_AlGaAs lambGap_AlGaAs] = Absorption(c,lambda,T) % Absorption(Al mole fraction in AlGaAs,wavelength [nm],temperature,)
%Absorption coefficient cauculus
%Results: [Absorption coefficient[m-1] Gap wavelength[nm]
 
Kb = 1.380648e-23; %Boltzmann Constant [J/K]
n = Dispersion_2(c,lambda,T);
 
lambda = (1e-9)*lambda; %[m]
 
 if c < 0.45
    %Gap Energy
    Egoq = (1.6e-19)*(1.519+(1.23*c)); %[J]
    aq = (1.6e-19)*(1e-4)*(4.9 + (0.7*c) + (3.7*(c^2)));
    tetaq = 250 + (5*c) + (260*(c^2));
    pq = 2.33+(1.7*c);
    eg_AlGaAs = (Egoq - (aq*tetaq*(nthroot((1+((2*T/tetaq)^pq)), pq)-1)/2));
    %Constants
    R_Bohr = 115e-10 - ((142e-10)*c) + ((61e-10)*(c^2)); % Exciton Bohr radius 
    Dieletrica = 13.18 - (3.12*c); % Dieletric constant
    ex_AlGaAs = 7.52e-22 + ((1.0912e-21)*c) + ((8.768e-22)*(c^2)); 
    lambGap_AlGaAs = (1.989e-25)/(eg_AlGaAs); %Conversion from energy[eV] to wavelength [m]
    % Absorption coefficient adjust parameters
    A = 3.1E-31 - (1.5825E-30)*c + ((1.6346E-29)*(c^2)) - ((7.575E-29)*(c^3)) + ((1.0542E-28)*(c^4));
    B = 1E-26 + ((1.6667E-26)*c) + ((1E-24)*(c^2)) - ((1.6667E-24)*(c^3)) ;
    C = - 2.2E-19 + ((3.425E-18)*c) - ((6.3708E-17)*(c^2)) + ((3.075E-16)*(c^3)) - ((4.2917E-16)*(c^4));
    E =  434.22528 + (382.70979*c) - (13419.67469*(c^2)) + (43571.66118*(c^3)) - (51704.54545*(c^4)) + (20795.62594*(c^5));
            if lambda<=lambGap_AlGaAs % Energia do foton > Energia do Gap
            % Teoria de Elliott
                zAlGaAs= sqrt(ex_AlGaAs/((1.989e-25/lambda) - eg_AlGaAs));
                a_AlGaAs = (((8*pi^2*Dieletrica)/((R_Bohr)^2*n*lambda))*(1/(1 - exp(-2*pi*zAlGaAs))))*(A*(1/lambda)^2 + B*(1/lambda) + C);
            else % Energia do foton < Energia do Gap
                n_Gap = Dispersion_2(c,(lambGap_AlGaAs)*1e9,T);
                D = ((8*pi^2*Dieletrica)/((R_Bohr)^2*n_Gap*lambGap_AlGaAs))*(A*(1/lambGap_AlGaAs)^2 + B*(1/lambGap_AlGaAs) + C);
                a_AlGaAs = D*exp(E*(((1.989e-25/lambda) - eg_AlGaAs)/(1.989e-25/lambda))*tanh((1.989e-25/lambda)/(2*Kb*300)));
            end
    else % 0.45 <= x < 1 (Indirect Gap)
        Egoq = (1.6e-19)*(1.519*(1-c) + 2.239*c - c*(1-c)*(-0.127+ 1.31*c));
        aq = (1.6e-19)*((4.9e-4)*(1-c) + (3.82e-4)*c);
        tetaq = (226)*(1-c) + (250)*c;
        pq =  (2.33)*(1-c) + (2.32)*c;
        eg_AlGaAs = (Egoq - (aq*tetaq*(nthroot((1+((2*T/tetaq)^pq)), pq)-1)/2));
        
        R_Bohr = 115e-10 - 142e-10*c + 61e-10*(c^2); % Exciton Bohr radius 
        Dieletrica = 13.18 - 3.12*c; % Dieletric constant
        ex_AlGaAs = 7.52e-22 + 1.0912e-21*c + 8.768e-22*(c^2); %Exciton binding energy
        lambGap_AlGaAs = (1.989e-25)/(eg_AlGaAs); %Conversion from energy[eV] to wavelength [m]
        % Absorption coefficient adjust parameters
        A = 1.3368E-31 - 3.7239E-31*c + 3.785E-31*(c^2) - 1.3679E-31*(c^3);
        B = - 2E-27 + 1E-26*c;
        C =  1.059E-18 - 6.0593E-18*c + 1.2235E-17*(c^2) - 1.0767E-17*(c^3) + 3.5E-18*(c^4);
        E =  434.22528 + 382.70979*c - 13419.67469*(c^2) + 43571.66118*(c^3) - 51704.54545*(c^4) + 20795.62594*(c^5);
        beta = 60E-20;
        if lambda<=lambGap_AlGaAs  % Energia do foton > Energia do Gap
            % Teoria de Elliott
            zAlGaAs = pi*sqrt((ex_AlGaAs)/((1.989e-25/lambda) - eg_AlGaAs));
            a_AlGaAs = ((2*beta)^2)*((8*pi^2*Dieletrica)/(3*(R_Bohr)^4*n*lambda))*((1 + (zAlGaAs^(-2)))/(1 - exp(-2*pi*zAlGaAs)));%*(A_k1*(1/lambda(b))^2 + B_k1*(1/lambda(b)) + C_k1);
        else
            n_Gap = Dispersion_2(c,(lambGap_AlGaAs)*1e9,T);
            D = (2*beta)^2*((8*pi^2*Dieletrica)/(3*(R_Bohr)^4*n_Gap*lambGap_AlGaAs))*(A*(1/lambGap_AlGaAs)^2 + B*(1/lambGap_AlGaAs) + C);
            a_AlGaAs = D*exp(E*(((1.989e-25/lambda) - eg_AlGaAs)/(1.989e-25/lambda))*tanh((1.989e-25/lambda)/(2*Kb*300)));
        end
 end
 lambGap_AlGaAs = lambGap_AlGaAs*(1e9); %in nanometers
end