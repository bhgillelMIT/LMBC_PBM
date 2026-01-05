function [DeltaH_rxn, H_products, H_reactants] = MP_Hrxn(Ts)

    %Function settings
    debug = false;
    
    if debug
        Ts = 300:100:2000;
    end


    %Define reference conditions
    Tref = 273.15 + 25; %K - reference state is 25 C
    pref = 101325; %Pa - reference state is 1 atm

%     %Specific heat of carbon
%     Cp_C = 10.68; %J/molK - Source: https://webbook.nist.gov/cgi/cbook.cgi?ID=C7440440&Units=SI&Mask=3 
%     
%     %Initial conditions
%     n_CH4 = 1;
%     X = 1; %Conversion
%     H0_CH4 = -74.8; %kJ/mol - Heat of Formation of Methane
%     Tf = 1200 + 273.15; %Final temperature of products
%     
%     n_C_f = X * n_CH4;
%     n_CH4_f = (1-X) * n_CH4;
%     n_H2_f = 2 * X * n_CH4;
%     
%     hi = H0_CH4;
%     Hf_C = n_C_f * Cp_C * (Tf - Tref);
%     Hf_H2 = n_H2_f * integral(@H2Cp, Tref, Tf);
%     Hf = Hf_C + Hf_H2
%     
%     
%     
%    
%     
%     if debug & false
%         
%          Ts = linspace(300, 2500);
%         Cps = H2Cp(Ts);
%         
%         figure()
%         
%         
%         plot(Ts, Cps);
%         
%         
%         
%         legend('Cp H_2');
%     end

%% From Heats of Formation and Specific Heats
    
    %Reference state
    T0 = 298.15;
    
    
    %Methane Properties
    CH4_DeltaH_form = -74870; %kJ/mol - Source: https://webbook.nist.gov/cgi/cbook.cgi?ID=C74828&Mask=1
    Cp_CH4 = @(T) CH4Cp(T)
    h_CH4 = @(T) CH4_DeltaH_form + integral(Cp_CH4, T0, T);
    
    %Hydrogen Properties 
    H2_DeltaH_form = 0;
    Cp_H2 = @(T) H2Cp(T)
    
    h_H2 = @(T) H2_DeltaH_form + integral(Cp_H2, T0, T);
    
    
    
    
    
    %Carbon Properties
    C_DeltaH_form = 0;
    C.A = 21.1751;
    C.B = -0.812428
    C.C = 0.448537
    C.D = -0.043256
    C.E = -0.013103   
    Cp_C = @(T) C.A + C.B.*(T./1000) + C.C .* (T./1000).^2 + C.D .* (T./1000).^3 + C.E./(T./1000).^2; 
    h_C = @(T) C_DeltaH_form + integral(Cp_C, T0, T);
    
    %Heat of Reaction 
    %DeltaH_rxn = 2 .* h_H2(T) + h_C(T) - h_CH4(T)
    
    for i = 1:length(Ts)
        
        T = Ts(i)
        
        
        H_products(i) = 2 .* h_H2(T) + h_C(T);
        H_reactants(i) = h_CH4(T);
        
        DeltaH_rxn(i) = H_products(i) - H_reactants(i);
        
        
        
        
    end
    
    
    
    

end



