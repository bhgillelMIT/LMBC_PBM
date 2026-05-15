%This function iterates through each cell and stores the downstream temperatures and conversions.

function params = IdentifyTempConvTracks(params, N_volumes_total)

    %Debug message
    fprintf('-- Logging downstream temperature and conversion.\n');

    %Settings
    sametol = 1E-6;

    %Allocate outputs
    params.Ts_above = cell(N_volumes_total,1); %Stores the temperatures above this bin
    params.Xs_above = cell(N_volumes_total,1); %Stores the conversions above this bin
    params.cell_active = cell(N_volumes_total,1); %Stores if the bin is trailing (false) or active (true)
    params.Ts_final = zeros(N_volumes_total,1);
    params.Xs_final = zeros(N_volumes_total, 1);

    %Only evaluate if heat is active
    if params.heat.active

        %Iterate through each volume 
        for iv = 1:N_volumes_total
    
            %Pull ind above  
            ind_above = iv; 
    
            %Iterate until you reach the final layer
            iz = params.zinds(iv);
            while ind_above > 0
    
    
                %Pull temperatures above
                T_above = params.Ts(ind_above);
                
    
                %Pull conversions above
                X_above = params.Xs(ind_above);
    
                %Determine if the bin is active here or not
                z_layer = params.zms(iz);
                T_largest = interp1(params.T_z.chars.largest.zs, params.T_z.chars.largest.Ts, z_layer);
                X_largest = interp1(params.T_z.chars.largest.zs, params.T_z.chars.largest.Xs, z_layer);
                X_match = (abs(X_above - X_largest)./X_largest < sametol);
                T_match = (abs(T_above - T_largest)./T_largest < sametol);
                active = ~(X_match & T_match); 
    
                %Store results
                params.Ts_above{iv}(end+1) = T_above;
                params.Xs_above{iv}(end+1) = X_above;
                params.cell_active{iv}(end+1) = active;

                %Get next ind above
                ind_above = params.inds_above{ind_above};
                iz = iz + 1;
    
            end

            %Store the final temperature and conversion
            if params.zinds(iv) == params.Nz
                params.Xs_final(iv) = params.Xs(iv);
                params.Ts_final(iv) = params.Ts(iv);
            else
                params.Xs_final(iv) = params.Xs_above{iv}(end);
                params.Ts_final(iv) = params.Ts_above{iv}(end);
            end

            
        end

        %Compile into arrays


    end


end