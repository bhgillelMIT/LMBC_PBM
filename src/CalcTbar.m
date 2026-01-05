function T_bar = CalcTbar(y_t, iz, im, params)

    

   

    if params.heat.active %Weighted average
        x = 1;      
    else %Assume uniform distribution of temperatures and conversions

        try

            %Pull appropriate temperature estimates
            Ts_cs = params.T_z.chars.Ts_cs{im};

            %Create a storage matrix
            T_mat = zeros(length(Ts_cs), params.chars.N_Xis);

            for it = 1:length(Ts_cs)
                Ts_c = Ts_cs{it};
                for ix = 1:params.chars.N_Xis
                    T_mat(it, ix) = Ts_c{ix}(iz);
        
                end
            end
    
            T_bar = mean(mean(T_mat));
        catch
            T_bar = params.T_liq;
        end

    end
    

end