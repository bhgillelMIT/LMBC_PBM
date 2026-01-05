function [Cps, hs] = CCp(Ts)

%% Define Functions 

    %Reference state
    T0 = 298.15;

    %Carbon Properties
    C_DeltaH_form = 0;
    C.A = 21.1751;
    C.B = -0.812428
    C.C = 0.448537
    C.D = -0.043256
    C.E = -0.013103   
    Cp_C = @(T) C.A + C.B.*(T./1000) + C.C .* (T./1000).^2 + C.D .* (T./1000).^3 + C.E./(T./1000).^2; 
    h_C = @(T) C_DeltaH_form + integral(Cp_C, T0, T);

%% Calculate all specific heats

    %Vector Inputs
    if exist('Ts', 'var')
        if length(Ts) > 1

            %Calculate specific heats
            Cps = Cp_C(Ts);
            
            for i = 1:length(Ts)
                T = Ts(i);
                hs(i) = h_C(T);
            end

        elseif length(Ts) == 1

            Cps = Cp_C(Ts);
            hs = h_C(Ts);
        else
            print('Ts does not exist')
        end
    
        
    else
        Cps = []
        
    end



end