function n = Dispersion(c,lambda,T)
    dispersion = 1; %Após análise dos resultados tops, percebi que o código 1 costuma ter resultados melhores.
    if(dispersion==1)
        %Sellmeier coefficients

        if c < 0.45 %Direct Gap
            A = 3.44 -2.8*c + 48*(c^2) -64*(c^3);
            B = 7.54 + 4.6*c - 64*(c^2) + 80*(c^3);
            C = 126424 - 187904*c + 1080735*(c^2) - 1185435*(c^3);
        else %Indirect Gap
            A = 14.7 - 22.7*c + 10*(c^2);
            B = -2.8 + 17*c - 8.8*(c^2);
            C = 313451 - 417013*c + 187003*(c^2);
        end

        n = sqrt(A + ((B*(lambda^2))/((lambda^2) - C))); %Sellmeier Equation


    elseif(dispersion==2)

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

end