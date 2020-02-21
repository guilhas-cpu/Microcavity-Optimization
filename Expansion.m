function eT = Expansion(s,e,c,Ta,T) % Expansion(thickness, Al mole fraction in AlGaAs, room temperature, desired temperature)

A = ((5.92e-7)*(exp(-6.88*c)))+3.18e-7;% Thermal expansion coefficient - Linear parameter ::: ((9.1e-7) - ((1.723e-6)*c));
B = -(((3.29e-12)*(exp(-5.54*c)))+1.75e-12);%Thermal expansion coefficient - Cubic parameter ::: (-5 + (8.733*c))*(1e-12);
  
%%%%%%%%%%%%%
%Correction of the thermal expansion coefficient due to refracive index
%%%%%%%%%%%%%
    if c == 0 %for GaAs only
        A = 7.55e-7;
        B = -3.396e-12;
    end
    
alfa = A*T + B*(T^3);

IdilTa = ((A/2)*(Ta^2))+((B/4)*(Ta^4));%Integral of thermal expansion coefficient at 300K
IdilT = ((A/2)*(T^2))+((B/4)*(T^4));%Integral of thermal expansion coefficient temperature dependent

Delta = (s/2)*e*(exp(IdilT) - exp(IdilTa)); %Expansion/Compression of the layer

eT = e + Delta; % final thisckness
end