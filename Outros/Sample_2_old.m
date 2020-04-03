function [xT xqwT xbT cavthickHalf cavthick xcavt] = Sample(s,ncs,nci,N,xqw,xb,c1,c2,c3,c4,cqw,lambdaR,T,teta,shift)
%sample(cavity type, # of QWs, QW thickness [nm], Barreir thickness [nm], Al mole fraction 1st
%layer, ri 2nd layer, ri 3rd layer, ri 4th layer, ressonant wavelength
%[nm], sample temperature [K])
global print;
global plota;

%Refractive index at 300K in the ressonant wavelength
n1R = Dispersion_2(c1,lambdaR,T); % Refractive index of the DBR 1st layer
n2R = Dispersion_2(c2,lambdaR,T); % Refractive index of the DBR 2nd layer
n3R = Dispersion_2(c3,lambdaR,T); % Refractive index of the cavity
n4R = Dispersion_2(c4,lambdaR,T); % Refractive index of the QWs, substrate and cap layer
[nqwR lambInGaAs EgInGaAs EgGaAs] = InGaAsDispersion(cqw,lambdaR,T);

Ta = 300; %Room temperature [K]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Building the sample
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%
%Defining the counters to build the DBR
%%%%
if n1R < n2R %If the 1st layer have a refractive index lower than the 2nd layer
    x(1) = 1; %Inserting a Caplayer with 1nm thickness [300K]
    xT(1) = Expansion(s,x(1),c4,Ta,T); %New thickness of caplayer
    beg = 2; %first counter for the 1st DBR with caplayer
    if ncs == 0
        fim = 2*nci+1; %counter for the lower DBR
    else
        if n3R >= n2R
            fim = 2*ncs+2; %final counter for the upper DBR
        else %n3R < n2R
            fim = 2*ncs+1; %final counter for the upper DBR
        end
    end
else %n1R >= n2R
    beg = 1; %remove the cap layer
    if ncs == 0
        fim = 2*nci; %counter for the lower DBR
    else
        if n3R < n2R
            fim = 2*ncs+1; %counter for the upper DBR
        else %n3R <= n2R
            fim = 2*ncs; %counter for the upper DBR
        end

    end
end

x1 = (lambdaR)/(4*n1R); %Thickness of 1st layer at 300K [nm] - from the air.
x1T = Expansion(s,x1,c1,Ta,T);%New thickness

x2 = (lambdaR)/(4*n2R); %Thickness of 2nd layer at 300K [nm].
x2T = Expansion(s,x2,c1,Ta,T);%New thickness

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%Building the Upper DBR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if n3R <= n2R
    for j = beg:fim
        if rem(j,2) == 0
            x(j)= x1; %Thickness of 1st layer (AlAs) at 300K [nm].
            xT(j) = x1T; %New thickness
        else
            x(j)= x2; % Thickness of 2nd layer (AlGaAs2) at 300K [nm].
            xT(j) = x2T; %New thickness
        end
    end
else %n3R > n2R
    for j = beg:fim
        if rem(j,2) == 0
            x(j)= x2; %Thickness of 1st layer (AlAs) at 300K [nm].
            xT(j) = x2T; %New thickness
        else
            x(j)= x1; % Thickness of 2nd layer (AlGaAs2) at 300K [nm].
            xT(j) = x1T; %New thickness
        end
    end
end

if nci == 0 || ncs == 0% For DBR only
    UpperDBRthick = sum(x);
    if(print)
        x
        disp('The upper DBR thickness at 300K is [nm]:');
        UpperDBRthick = UpperDBRthick
    end
    xqwT = 0;
    xbT = 0;
    cavthickHalf = 0;
    cavthick = 0;
    xcavt = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Refractive index pattern
%%%%%%%%%%%%%%%%%%%%%%%%%%%

Air = 100; %Air gap [nm]
Sub = 500; %Subtrate thickness [nm]
np = 35000; %number of points
t = sum(xT);

z = linspace(0,t+Air+Sub,np); % Growth axes [nm]

for b = 1:np 
    if z(b) <= Air %The air gap layer
       ri(b) = 1; %Refractive index of the air
       fim = b;
       esp = z(b);
    else
        break
    end
end

%nl = n1R; % Refractive index of the DBR 1st layer;

for d = 1:length(xT)% Scaning each layer
    g = esp + xT(d);
    
    if xT(d) == xT(1) %Cap layer
        for b = (fim+1):np 
            if z(b) <= g 
                ri(b) = n4R; %GaAs
                fim = b;
                esp = z(b);
            else
                break
            end
        end
    end
    
    if xT(d) == x1T %1st layer
        for b = (fim+1):np %Scan the sample
            if z(b) <= g
                ri(b) = n1R; %Cavity layer refractive index
                esp = z(b);
                fim = b;
            else
                break
            end
        end
    end
    
    if xT(d) == x2T %2nd layer
        for b = (fim+1):np %Scan the sample
            if z(b) <= g
                ri(b) = n2R; %Cavity layer refractive index
                esp = z(b);
                fim = b;
            else
                break
            end
        end
    end
    
    if xT(d) == cavthick %Cavity layer 
        for b = (fim+1):np %Scan the sample
            if z(b) <= g
                ri(b) = n3R; %Cavity layer refractive index
                esp = z(b);
                fim = b;
            else
                break
            end
        end
    end 
    
    if xT(d) == cavthickHalf %Cavity layer cavthickHalf
        for b = (fim+1):np %Scan the sample
            if z(b) <= g
                ri(b) = n3R; %Cavity layer refractive index
                esp = z(b);
                fim = b;
            else
                break
            end
        end
    end 
    
    if xT(d) == xqwT %QW layer
        for b = (fim+1):np %Scan the sample
            if z(b) <= g
                ri(b) = n4R; %QW refractive index
                esp = z(b);
                fim = b;
            else
                break
            end
        end
    end
    
    if xT(d) == xbT %Barreir layer
        for b = (fim+1):np %Scan the sample
            if z(b) <= g
                ri(b) = n3R; %Cavity layer refractive index
                esp = z(b);
                fim = b;
            else
                break
            end
        end
    end
end

for b = (fim+1):np %Substrate
    if z(b) <= t+Air+Sub
        ri(b) = n4R;
    else
        break
    end
end

if(plota)
    figure()
    plot(z,ri)
    xlabel('z(nm)')
    ylabel('Refractive Index')
    title('Refractive Index Pattern','FontSize',10)
end
return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Building the cavity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if s == 1 %lambda/2 cavity
    if n3R <= n2R
        if N == 0
            if(print)
                disp('Passive cavity!');
            end
            n = length(x);
            x(n+1) = s*lambdaR/(2*n3R) + shift;%Cavity
            xT(n+1) = Expansion(s,x(n+1),c3,Ta,T);%New thickness
            cavthick = xT(n+1);
            xcavt = cavthick;
            xb = 0;
            xbT = 0;
            xqwT = 0;
            cavthickHalf = 0;
        else % N>0
            if N == 1
                xb = 0;
                xbT = 0;
            end
            n = length(x);
            x(n+1) = (lambdaR/(2*n3R) - (N*xqw + (N-1)*xb) + shift)/2 ;% 1st part of the cavity
            xT(n+1) = Expansion(s,x(n+1),c3,Ta,T);%New thickness
            cavthick = xT(n+1);
            cavthickHalf = 0;
            n = length(x);
            for j = 1:((2*N)-1) % Building the actve medium
                if rem(j,2) == 0 % Barreir
                    x(n+j) = xb; 
                    xT(n+j)= Expansion(s,x(n+j),c3,Ta,T);%New thickness
                    xbT = xT(n+j);
                else %QW
                    x(n+j) = xqw; %QWs GaAs thickness [nm]
                    xT(n+j) = Expansion(s,x(n+j),c4,Ta,T);%New thickness
                    xqwT = xT(n+j);
                end
            end
            n = length(x);
            x(n+1) = (lambdaR/(2*n3R) - (N*xqw + (N-1)*xb) + shift)/2;% 2st part of the cavity 
            xT(n+1) = Expansion(s,x(n+1),c3,Ta,T);%New thickness
            xcavt = (2*xT(n+1)) + N*xqwT + (N-1)*xbT;
        end
    else %n3R > n2R
         if N == 0
            if(print)
                disp('Passive cavity!');
            end
            n = length(x);
            x(n+1) = s*lambdaR/(2*n3R) + shift;%Cavity
            xT(n+1) = Expansion(s,x(n+1),c3,Ta,T);%New thickness
            cavthick = xT(n+1);
            xb = 0;
            xbT = 0;
            xqw = 0;
            xqwT = 0;
            xcavt = cavthick;
            cavthickHalf = 0;
         else
            if N == 1 
                n = length(x);
                x(n+1) = xb; %Barreir thickness [nm]
                xT(n+1) = Expansion(s,x(n+1),c3,Ta,T);%New thickness
                
                x(n+2) = xqw;
                xT(n+2) = Expansion(s,x(n+2),c4,Ta,T);%New thickness
                
                x(n+3) =(lambdaR/(2*n3R)) - 2*(xqw+xb) + shift;% s=1 cavity
                xT(n+3) = Expansion(s,x(n+3),c3,Ta,T);%New thickness
                cavthick = xT(n+3);
                cavthickHalf = 0;
                
                x(n+4) = xqw; %QWs GaAs thickness [nm]
                xT(n+4) = Expansion(s,x(n+4),c4,Ta,T);%New thickness
                xqwT = xT(n+4);
                
                x(n+5) = xb; %Barreir thickness [nm]
                xT(n+5)= Expansion(s,x(n+5),c3,Ta,T);%New thickness 
                xbT = xT(n+5);
                
                xcavt = xT(n+1) + xT(n+2) + xT(n+3) + xT(n+4) + xT(n+5);
            else %N>1
                n = length(x);
                for j = 1:(2*N) % Building the actve medium on the edge of the DBR
                    if rem(j,2) == 0 %QW
                        x(n+j) = xqw; %QWs GaAs thickness [nm]
                        xT(n+j) = Expansion(s,x(n+j),c4,Ta,T);%New thickness
                        xqwT = xT(n+j);
                    else 
                        x(n+j) = xb; % Barreir
                        xT(n+j) = Expansion(s,x(n+j),c3,Ta,T);%New thickness
                        xbT = xT(n+j);
                    end
                end
                n = length(x);
                x(n+1) = lambdaR/(2*n3R) - 2*N*(xqw + xb) + shift;% s=1 cavity
                xT(n+1) = Expansion(s,x(n+1),c3,Ta,T);%New thickness
                cavthick = xT(n+1);
                n = length(x);
                for j = 1:(2*N) % Building the actve medium on the edge of the DBR
                    if rem(j,2) == 0 
                        x(n+j) = xb; % Barreir
                        xT(n+j)= Expansion(s,x(n+j),c3,Ta,T);%New thickness
                        xbT = xT(n+j);
                    else %QW
                        x(n+j) = xqw; %QWs GaAs thickness [nm]
                        xT(n+j)= Expansion(s,x(n+j),c4,Ta,T);%New thickness
                        xqwT = xT(n+j);
                    end
                end
                xcavt = (2*N*(xqwT + xbT)) + cavthick;
                cavthickHalf = 0;
            end
        end
    end
else %s>1
    if N == 0
        if(print)
            disp('Passive cavity!');
        end
        n = length(x);
        x(n+1) = s*lambdaR/(2*n3R) + shift/s;%Cavity
        xT(n+1) = Expansion(s,x(n+1),c3,Ta,T);%New thickness
        cavthick = xT(n+1);
        cavthickHalf = 0;
        xbT = 0;
        xqwT = 0;
        xcavt = cavthick;
    else % Cavity with QWs
        if n1R > n3R 
            if N == 1
                xbT = 0;
                n = length(x);
                x(n+1) = ((lambdaR/(2*n3R)) - xqw + shift/s)/2;% First part of cavity
                xT(n+1) = Expansion(s,x(n+1),c3,Ta,T);%New thickness
                cavthickHalf = xT(n+1);
                n = length(x);
                for j = 1:(s-1)
                    x(n+j) = xqw; %SQW GaAs thickness [nm]
                    xT(n+j)= Expansion(s,x(n+j),c4,Ta,T);%New thickness
                    xqwT = xT(n+j);
                    
                    n = length(x);
                    x(n+1) = (lambdaR/(2*n3R)) - xqw + shift/s; %Second part of cavity 
                    xT(n+1) = Expansion(s,x(n+1),c3,Ta,T);%New thickness
                    cavthick = xT(n+1);
                    n = length(x)-1;
                end
                n = length(x);
                x(n+1) = xqw; % Last SQW GaAs thickness [nm]
                xT(n+1) = Expansion(s,x(n+1),c4,Ta,T);%New thickness
                xqwT = xT(n+1);
                
                x(n+2) = ((lambdaR/(2*n3R)) - xqw + shift/s)/2;% Last part of cavity
                xT(n+2) = Expansion(s,x(n+2),c3,Ta,T); %New thickness
                
                xcavt = (2*cavthickHalf) + ((s-1)*(xqwT + cavthick)) + xqwT;
            else %N>1 %aqui!!!!
                n = length(x);
                x(n+1) = ((lambdaR/(2*n3R)) - ((s-1)*(N*xqw + (N-1)*xb)) + (shift/s))/2;% First part of cavity
                xT(n+1) = Expansion(s,x(n+1),c3,Ta,T); %New thickness
                cavthickHalf = xT(n+1);
                n = length(x);
                for j = 1:(s-1) %Active medium
                    for j = 1:((2*N)-1) % Building the actve medium
                        if rem(j,2) == 0 % Barreir
                            x(n+j) = xb;
                            xT(n+j) = Expansion(s,x(n+j),c3,Ta,T);%New thickness
                            xbT = xT(n+j);
                        else %QW
                            x(n+j) = xqw; %QWs GaAs thickness [nm]
                            xT(n+j) = Expansion(s,x(n+j),c4,Ta,T);%New thickness
                            xqwT = xT(n+j);
                        end
                    end
                    n = length(x);
                    x(n+1) = (lambdaR/(2*n3R)) - ((s-1)*(N*xqw + (N-1)*xb)) + (shift/s);% 2nd part of cavity
                    xT(n+1) = Expansion(s,x(n+1),c3,Ta,T);%New thickness
                    cavthick = xT(n+1);
                    n = length(x);
                end
                n = length(x);
                for j = 1:((2*N)-1) % Building the actve medium
                    if rem(j,2) == 0 % Barreir
                        x(n+j) = xb;
                        xT(n+j) = Expansion(s,x(n+j),c3,Ta,T);%New thickness
                        xbT = xT(n+j);
                    else %QW
                        x(n+j) = xqw; %QWs GaAs thickness [nm]
                        xT(n+j) = Expansion(s,x(n+j),c4,Ta,T);%New thickness
                        xqwT = xT(n+j);
                    end
                end
                n = length(x);
                x(n+1) = ((lambdaR/(2*n3R)) - ((s-1)*(N*xqw + (N-1)*xb)) + (shift/s))/2;% First part of cavity
                xT(n+1) = Expansion(s,x(n+1),c3,Ta,T); %New thickness
                
                xcavt = (2*cavthickHalf) + s*(N*xqwT + (N-1)*xbT) + (s-1)*(cavthick);
            end 
        else %n1R < n3R
            if N == 1 
                n = length(x);
                for j = 1:(s-1)
                    x(n+j) = xb; %Barreir thickness [nm]
                    xT(n+j) = Expansion(s,x(n+j),c3,Ta,T);%New thickness
                    xbT = xT(n+j);
                    
                    n = length(x);
                    
                    x(n+1) = (lambdaR/(2*n3R)) - xqw + shift/s; %Cavity Part
                    xT(n+1)= Expansion(s,x(n+1),c3,Ta,T);%New thickness
                    cavthick = xT(n+1);
                    n = length(x)-1;
                end
                x(n+1) = ((lambdaR/(2*n3R)) - xqw + shift/s)/2;% %Cavity Part
                xT(n+1) = Expansion(s,x(n+1),c3,Ta,T);%New thickness
                cavthickHalf = xT(n+1);
                
                n = length(x);
                x(n+1) = xqw; % Last SQW GaAs thickness [nm]
                xT(n+1) = Expansion(s,x(n+1),c4,Ta,T);%New thickness
                xqwT = xT(n+1);
                
                x(n+2) = ((lambdaR/(2*n3R)) - xqw + shift/s)/2;% %Last cavity Part
                xT(n+2)= Expansion(s,x(n+2),c3,Ta,T);%New thickness
                
                xcavt = (s-1)*(xbT + cavthick) + 2*cavthickHalf + xqwT
            else %N>1    
                n = length(x);
                for j = 1:(2*N) %Building the actve medium on the edge of the DBR
                    if rem(j,2) == 0 %QW
                        x(n+j) = xqw; %QWs GaAs thickness [nm]
                        xT(n+j)= Expansion(s,x(n+j),c4,Ta,T);%New thickness
                        xqwT = xT(n+j);
                    else 
                        x(n+j) = xb; % Barreir
                        xT(n+j) = Expansion(s,x(n+j),c3,Ta,T);%New thickness
                        xbT = xT(n+j);
                    end
                end
                for b = 1:(s-1);%Building the actve medium inside the cavity
                    n = length(x);
                    x(n+1) = (lambdaR/(2*n3R)) - (((s-1)/s)*(N*xqw + (N-1)*xb)) - (2*N*(xqw+xb)/s + shift/s);%First cavity part
                    xT(n+1) = Expansion(s,x(n+1),c3,Ta,T);%New thickness
                    cavthick = xT(n+1);
                    cavthickHalf = 0;
                    n = length(x);
                    for j = 1:((2*N)-1) % Building the actve medium inside the cavity
                        if rem(j,2) == 0 % Barreir
                            x(n+j) = xb;
                            xT(n+j) = Expansion(s,x(n+j),c3,Ta,T);%New thickness
                            xbT = xT(n+j);
                        else %QW
                            x(n+j) = xqw; %QWs GaAs thickness [nm]
                            xT(n+j) = Expansion(s,x(n+j),c4,Ta,T);%New thickness
                            xqwT = xT(n+j);
                        end
                    end
                end
                n = length(x);
                x(n+1) = (lambdaR/(2*n3R)) - (((s-1)/s)*(N*xqw + (N-1)*xb)) - (2*N*(xqw+xb)/s + shift/s);%Last cavity part
                xT(n+1) = Expansion(s,x(n+1),c3,Ta,T);%New thickness
                n = length(x);
                for j = 1:(2*N) %Building the actve medium on the edge of the DBR
                    if rem(j,2) == 0 %QW
                        x(n+j) = xb; % Barreir
                        xT(n+j) = Expansion(s,x(n+j),c3,Ta,T);%New thickness
                        xbT = xT(n+j);
                    else 
                        x(n+j) = xqw; %QWs GaAs thickness [nm]
                        xT(n+j) = Expansion(s,x(n+j),c4,Ta,T);%New thickness
                        xqwT = xT(n+j);
                    end
                end
            xcavt = 2*N*(xqwT + xbT) + (s-1)*(N*xqwT + (N-1)*xbT) + s*cavthick;
            end
        end 
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Bottom DBR
%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = length(x); 
   %Constructing the botom DBR
    if n3R <= n2R
        %disp('OK!');
        xe = x1;
        xeT = x1T;
        for j = 1 : (2*nci)+1%Extra layer
            if rem(j,2) == 0      
                x(n+j)= xe; %Thickness of 1st layer (AlAs) at 300K [nm].
                xT(n+j) = xeT; %New thickness
                if xe == x1 %Exchanging the refractive index
                    xe = x2;
                    xeT = x2T;
                else
                    xe = x1;
                    xeT = x1T;
                end
            else
                x(n+j)=xe; % Thickness of 1st layer at 300K [nm].
                xT(n+j) = xeT; %New thickness
                if xe == x1 %Exchanging the refractive index
                    xe = x2;
                    xeT = x2T;
                else
                    xe = x1;
                    xeT = x1T;
                end
            end
        end
    else %n3R > n1R
        xe = x2;
        xeT = x2T;
        for j = 1 : (2*nci)+1%Extra layer
            if rem(j,2) == 0      
                x(n+j)= xe; %Thickness of 1st layer (AlAs) at 300K [nm].
                xT(n+j) = xeT; %New thickness
                if xe == x1 %Exchanging the refractive index
                    xe = x2;
                    xeT = x2T;
                else
                    xe = x1;
                    xeT = x1T;
                end
            else
                x(n+j)=xe; % Thickness of 1st layer (AlGaAs2) at 300K [nm].
                xT(n+j) = xeT; %New thickness
                if xe == x1 %Exchanging the refractive index
                    xe = x2;
                    xeT = x2T;
                else
                    xe = x1;
                    xeT = x1T;
                end
            end
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Refractive index pattern
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Air = 100; %Air gap [nm]
Sub = 500; %Subtrate thickness [nm]
np = 35000; %number of points
t = sum(xT);

z = linspace(0,t+Air+Sub,np); % Growth axes [nm]

for b = 1:np 
    if z(b) <= Air %The air gap layer
       ri(b) = 1; %Refractive index of the air
       fim = b;
       esp = z(b);
    else
        break
    end
end

%nl = n1R; % Refractive index of the DBR 1st layer;

for d = 1:length(xT)% Scaning each layer
    g = esp + xT(d);
    
    if xT(d) == xT(1) %Cap layer
        for b = (fim+1):np 
            if z(b) <= g 
                ri(b) = n4R; %GaAs
                fim = b;
                esp = z(b);
            else
                break
            end
        end
    end
    
    if xT(d) == x1T %1st layer
        for b = (fim+1):np %Scan the sample
            if z(b) <= g
                ri(b) = n1R; %Cavity layer refractive index
                esp = z(b);
                fim = b;
            else
                break
            end
        end
    end
    
    if xT(d) == x2T %2nd layer
        for b = (fim+1):np %Scan the sample
            if z(b) <= g
                ri(b) = n2R; %Cavity layer refractive index
                esp = z(b);
                fim = b;
            else
                break
            end
        end
    end
    
    if xT(d) == cavthick %Cavity layer 
        for b = (fim+1):np %Scan the sample
            if z(b) <= g
                ri(b) = n3R; %Cavity layer refractive index
                esp = z(b);
                fim = b;
            else
                break
            end
        end
    end 
    
    if xT(d) == cavthickHalf %Cavity layer cavthickHalf
        for b = (fim+1):np %Scan the sample
            if z(b) <= g
                ri(b) = n3R; %Cavity layer refractive index
                esp = z(b);
                fim = b;
            else
                break
            end
        end
    end 
    
    if xT(d) == xqwT %QW layer
        for b = (fim+1):np %Scan the sample
            if z(b) <= g
                ri(b) = nqwR; %QW refractive index
                esp = z(b);
                fim = b;
            else
                break
            end
        end
    end
    
    if xT(d) == xbT %Barreir layer
        for b = (fim+1):np %Scan the sample
            if z(b) <= g
                ri(b) = n3R; %Cavity layer refractive index
                esp = z(b);
                fim = b;
            else
                break
            end
        end
    end
end

for b = (fim+1):np %Substrate
    if z(b) <= t+Air+Sub
        ri(b) = n4R;
    else
        break
    end
end

if(plota)
    figure()
    plot(z,ri)
    xlabel('z(nm)')
    ylabel('Refractive Index')
    title('Refractive Index Pattern','FontSize',10)
end

xT = round(xT,1);
xqwT = round(xqwT,1);
xbT = round(xbT,1);
cavthickHalf = round(cavthickHalf,1);
cavthick = round(cavthick,1);
xcavt = round(xcavt,1);
    x = round(x,1);
    xT = round(xT,1);
    if(print)
        x = x; %disabled 
        disp('Microcavity total thickness at 300K [nm]:');
        sum(x)
        xT = xT; %disabled
        disp('Microcavity total thickness at defined temperature [nm]:');
        sum(xT) 
    end
end
