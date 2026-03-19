function X_bar = CalcXbar(y_t, iz, im, params)

    

    if params.react.active %Weighted average
        x = 1;      
    else %Assume uniform distribution of temperatures and conversions

        try

            %Pull appropriate temperature estimates
            Xs_cs = params.T_z.chars.Xs_cs{im};

            %Pull trailing inds
            inds_trail = params.T_z.chars.itrails{im};

            %Identify active bins
            if iz == 1
                active_inds = params.T_z.chars.merged{im}.bins_in(:,1);
            else
                active_inds = find(params.T_z.chars.merged{im}.Tcz(:,iz) > 0); 
            end
            active_inds = active_inds(active_inds > max(inds_trail));
            active_inds = active_inds - 1;
        
            %Create a storage matrix
            X_mat = zeros(length(active_inds), params.chars.N_Xis);

            for it = 1:length(active_inds) %1:length(Xs_cs)
                ind = active_inds(it);
                Xs_c = Xs_cs{ind};
                for ix = 1:params.chars.N_Xis
                    X_mat(it, ix) = Xs_c{ix}(iz);
        
                end
            end
    
            %Calculate overall mean
            X_bar = mean(X_mat(:,1)); %mean(mean(X_mat));

        catch
            X_bar = 0;
        end
    end


end