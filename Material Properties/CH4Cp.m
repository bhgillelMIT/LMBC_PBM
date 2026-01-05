function Cps = CH4Cp(Ts)

%% Function Setup

    debug = false;
    

%% Define functions 

    %Hydrogen specific heat func - Source: https://webbook.nist.gov/cgi/cbook.cgi?ID=C74828&Mask=1
    CH4.A1300 = -0.703029;
    CH4.B1300 = 108.4773;
    CH4.C1300 = -42.52157;
    CH4.D1300 = 5.862788;
    CH4.E1300 = 0.678565;
    CH4.F1300 = -76.84376;
    CH4.G1300 = 158.7163;
    CH4.H1300 = -74.87310;
    CH4.Cp1300 = @(T) CH4.A1300 + CH4.B1300 .* (T./1000) + CH4.C1300 .* (T./1000).^2 ...
                     + CH4.D1300 .* (T./1000).^3 + CH4.E1300./((T./1000).^2); 
                   
    
    CH4.A6000 = 85.81217;
    CH4.B6000 = 11.26467;
    CH4.C6000 = -2.114146;
    CH4.D6000 = 0.138190;
    CH4.E6000 = -26.42221;
    CH4.F6000 = -153.5327;
    CH4.G6000 = 224.4143;
    CH4.H6000 = -74.87310;
    CH4.Cp6000 = @(T) CH4.A6000 + CH4.B6000 .* (T./1000) + CH4.C6000 .* (T./1000).^2 ...
                     + CH4.D6000 .* (T./1000).^3 + CH4.E6000./((T./1000).^2); 
    
    
    
%% Debug plots

    if debug
        Ts1000 = linspace(300,1300);
        Cps1300 = CH4.Cp1300(Ts1000);
        Ts2500 = linspace(1300,2500);
        Cps6000 = CH4.Cp6000(Ts2500);
        figure()
        plot(Ts1000, Cps1300); hold on;
        plot(Ts2500, Cps6000);
        title('CH_{4} Cp');
        xlabel('Temperature (K)');
        ylabel('Specific Heat (J/molK)')
    end

%% Calculate all specific heats

    %Vector Inputs
    if length(Ts) > 1

        %Identify which inputs fall into which temp range
        sub1300s = Ts <= 1300;
        sub6000s = Ts >= 1300 & Ts <= 2500;

        %Allocate output vector
        Cps = zeros(size(Ts));

        %Calculate specific heats
        Cps(sub1300s) = CH4.Cp1300(Ts(sub1300s));
        Cps(sub6000s) = CH4.Cp6000(Ts(sub6000s));
    
    else
        
        if Ts <= 1300
            Cps = CH4.Cp1300(Ts);
        elseif Ts <=6000
            Cps = CH4.Cp6000(Ts);
        else
            error('Invalid Temperature Input (exceeds 6000 K)')
        end
    end
    



end