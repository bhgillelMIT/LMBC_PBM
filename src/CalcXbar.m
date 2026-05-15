function [X_bar, Xs_m, Ns_m] = CalcXbar(y, y_t, iz, im, params)

    %Pull overall cell value
    ind_m = find(params.zinds_m(:) == iz & params.xinds_m(:) == im);
    xinds_zm = find(params.xinds == im & params.zinds == iz);
    
    

    if params.react.active %Weighted average

        %Pull values
        ys_zm = y(xinds_zm);
        Xs_m = params.Xs(xinds_zm);
        %Ts_m = params.Ts(xinds_zm);

        %Calculate y fraction
        y_frac = ys_zm/(sum(ys_zm)+1E-32);

        %Calcualte weighted average
        if any(y_frac > 0)
            X_bar = sum(y_frac(:) .* Xs_m(:));
        else
            X_bar = mean(Xs_m);
        end

        %Define other outputs
        Ns_m = ys_zm;


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

            %Output the distribution
            Xs_m = 0; %TO BE COMPLETED
            Ns_m = 0;

        catch
            X_bar = 0; Xs_m = 0; Ns_m = 0;
        end
    end


end