function Cps = H2Cp(Ts)

%% Function Setup

    debug = false;
    

%% Define functions 

    %Hydrogen specific heat func - Source: https://webbook.nist.gov/cgi/cbook.cgi?ID=C1333740&Mask=1&Type=JANAFG&Table=on
    H2.A1000 = 33.066178;
    H2.B1000 = -11.363417;
    H2.C1000 = 11.432816;
    H2.D1000 = -2.772874;
    H2.E1000 = -0.158558;
    H2.F1000 = -9.980797;
    H2.G1000 = 172.707974;
    H2.H1000 = 0;
    H2.Cp1000 = @(T) H2.A1000 + H2.B1000 .* (T./1000) + H2.C1000 .* (T./1000).^2 ...
                     + H2.D1000 .* (T./1000).^3 + H2.E1000./((T./1000).^2); 
                   
    
    H2.A2500 = 18.5630838;
    H2.B2500 = 12.257357;
    H2.C2500 = -2.859786;
    H2.D2500 = 0.268238;
    H2.E2500 = 1.977990;
    H2.F2500 = -1.147438;
    H2.G2500 = 156.288133;
    H2.H2500 = 0;
    H2.Cp2500 = @(T) H2.A2500 + H2.B2500 .* (T./1000) + H2.C2500 .* (T./1000).^2 ...
                     + H2.D2500 .* (T./1000).^3 + H2.E2500./((T./1000).^2); 
    
    
    
%% Debug plots

    if debug
        Ts1000 = linspace(300,1000);
        Cps1000 = H2.Cp1000(Ts1000);
        Ts2500 = linspace(1000,2500);
        Cps2500 = H2.Cp2500(Ts2500);
        figure()
        plot(Ts1000, Cps1000); hold on;
        plot(Ts2500, Cps2500);
        title('H_{2} C_p');
        xlabel('Temperature (K)');
        ylabel('Specific Heat (J/molK)')
    end

%% Calculate all specific heats

    %Vector Inputs
    if exist('Ts', 'var')
        if length(Ts) > 1

            %Identify which inputs fall into which temp range
            sub1000s = Ts <= 1000;
            sub2500s = Ts >= 1000 & Ts <= 2500;

            %Allocate output vector
            Cps = zeros(size(Ts));

            %Calculate specific heats
            Cps(sub1000s) = H2.Cp1000(Ts(sub1000s));
            Cps(sub2500s) = H2.Cp2500(Ts(sub2500s));


        elseif length(Ts) == 1

            if Ts <= 1000
                Cps = H2.Cp1000(Ts);
            elseif Ts <=2500
                Cps = H2.Cp2500(Ts);
            else
                error('Invalid Temperature Input (exceeds 2500 K)')
            end

        else
            print('Ts does not exist')
        end
    
        
    else
        Cps = []
        
    end


end