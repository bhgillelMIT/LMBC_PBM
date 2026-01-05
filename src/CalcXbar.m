function X_bar = CalcXbar(y_t, iz, im, params)

    

    if params.react.active %Weighted average
        x = 1;      
    else %Assume uniform distribution of temperatures and conversions

        try

            %Pull appropriate temperature estimates
            Xs_cs = params.T_z.chars.Xs_cs{im};
        
            %Create a storage matrix
            X_mat = zeros(length(Xs_cs), params.chars.N_Xis);

            for it = 1:length(Xs_cs)
                Xs_c = Xs_cs{it};
                for ix = 1:params.chars.N_Xis
                    X_mat(it, ix) = Xs_c{ix}(iz);
        
                end
            end
    
            X_bar = mean(mean(X_mat));

        catch
            X_bar = 0;
        end
    end


end