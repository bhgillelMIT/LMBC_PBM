%This function should be run in PBM_v3 before PBM_solver is run. Can then
%remove the calculation of Ts from the simulation since they're known. Only
%the numeric densities (y) need to be updated each iteration.

function params = CalcCharTemps(y, params)

    if params.heat.active


        %Iterate through z and m
        for iz = 1:params.Nz
            for im = 1:params.Nms
    
                %Find corresponding indices
                relinds = find(params.xinds == im & params.zinds == iz);

                %Pull numeric densities
                relys = y(relinds); 

                %Handle initial empty case
                if all(relys == 0)
                    relys = ones(size(relys));
                end

                %Calculate mean temperature
                relys_sum = sum(relys);
                if relys_sum > 0 %Calculate weighted mean (not mixing cup)
                    yfrac = relys/relys_sum;
                    if params.react.active
                        [Ts,Xs] = ExtractTempsConvs(iz, im, params);
                        T_bar = sum(yfrac(:) .* Ts(:));
                        x = 1; %To be implemented
                    else
                        Ts = ExtractTemps(iz, im, 1, params);
                        T_bar = sum(yfrac(:) .* Ts(:));
                    end
                    
                else %Just use an average of the temperature points temporarily

                    if params.react.active
                        x = 1; %To be implemented
                    else
                        Ts = ExtractTemps(iz, im, 1, params);
                        T_bar = mean(Ts);
                    end
                end

                if T_bar < 0
                    x = 1;
                end

                %Log result
                params.Xs(relinds) = Xs;
                params.Ts(relinds) = Ts;
                params.T_mu(iz,im) = T_bar;
    
    
            end
        end


    end



end



function [Tsm] = ExtractTemps(iz, im, ix, params)

    %Pull the initial flux 
    init_flux = params.N_dot_o(im);

    



    %Case 1 - Initial layer
    if iz == 1

        %Iterate through and extract temperatures
        NT = length(params.T_z.chars.Ts_cs{im});
        Tsm = zeros(1, NT);
        bins_active = 1:NT;
        for it = 1:NT
            if it == NT
                %Calculate average between final pivot and liquid temperature
                Tsm(it) = (params.T_z.chars.Ts_cs{im}{it}{ix}(iz) + params.Tsz(iz))./2; %Use liquid temperature as upper bound
            else
                Tsm(it) = (params.T_z.chars.Ts_cs{im}{it}{ix}(iz) + params.T_z.chars.Ts_cs{im}{it+1}{ix}(iz))./2;
            end
        end

        %Calculate normal numeric densities
        %init_fracs = params.T_z.chars.fracs_in{im};
        %init_fluxes = init_flux .* init_fracs;


    %Case 2 - Further layers 
    else
        merged = params.T_z.chars.merged{im}{1};
        bins_in = find(merged.Tcz(:,iz-1) > -1);
        bins_out = find(merged.Tcz(:,iz) > -1);
        Tsb = merged.Tcz(bins_out, iz-1);
        Tsm = (Tsb(1:end-1) + Tsb(2:end))./2;
        
    end


    % bins_active = [bins_active, NT+1]; %Account for liquid temperature bound
    %     Tbs = merged.Tcz(bins_active,iz);
    %     Tcs = 1;
    %     Ts(it) = paramas.T_z.chars.merged
    % Ts = 1;


end



function [Tsm, Xsm] = ExtractTempsConvs(iz, im, params)

    NX = params.chars.N_Xis;

    %Calculate value for the largest
    T_largest = interp1(params.T_z.chars.largest.zs, params.T_z.chars.largest.Ts, params.zms(iz));
    X_largest = interp1(params.T_z.chars.largest.zs, params.T_z.chars.largest.Xs, params.zms(iz));


    %Case 1 - Initial layer
    if iz == 1

        %Iterate through and extract temperatures
        NT = length(params.T_z.chars.Ts_cs{im});
        
        Tsm = zeros(1, NT*NX);
        bins_active = 1:NT;
        iters = 1;
        for it = 1:NT
            for ix = 1:NX

                %Calculate mean temperature
                if it == NT
                    %Calculate average between final pivot and liquid temperature
                    Tsm(iters) = (params.T_z.chars.Ts_cs{im}{it}{ix}(iz) + params.Tsz(iz))./2; %Use liquid temperature as upper bound
                else
                    Tsm(iters) = (params.T_z.chars.Ts_cs{im}{it}{ix}(iz) + params.T_z.chars.Ts_cs{im}{it+1}{ix}(iz))./2;
                end
                

                %Calculate mean conversion
                if ix == NX
                    Xsm(iters) = (params.T_z.chars.Xs_cs{im}{it}{ix}(iz) + 1)./2;
                else
                    Xsm(iters) = (params.T_z.chars.Xs_cs{im}{it}{ix}(iz) + params.T_z.chars.Xs_cs{im}{it}{ix+1}(iz))./2;
                end

                %Update iteration counter
                iters = iters + 1;

            end

            

        end

        %Calculate normal numeric densities
        %init_fracs = params.T_z.chars.fracs_in{im};
        %init_fluxes = init_flux .* init_fracs;


    %Case 2 - Further layers 
    else
        
        layer_inds = (params.zinds == iz & params.xinds == im);
        Tinds = params.Tinds(layer_inds);
        Xinds = params.Xinds(layer_inds);

        Tsm = zeros(size(Tinds));
        Xsm = zeros(size(Tinds));

        
        for ix = 1:params.chars.N_Xis


            %Identify bins entering and exiting
            merged = params.T_z.chars.merged{im}{ix};
            Tcz = merged.Tcz;
            Tbz = merged.Tbz;
            Xcz = merged.Xcz;
            Xbz = merged.Xbz;
            bins_in = find(merged.Tcz(:,iz-1) > -1);
            bins_out = find(merged.Tcz(:,iz) > -1);
    
            NT = length(bins_out);
            %Tsm = zeros(1, NT*NX);
            %Iterate through conversions
            Tsb = merged.Tcz(bins_out, iz);
            Xsb = merged.Xcz(bins_out, iz);
            for it = 1:(NT-1)

                %Calculate index
                ind = find(Tinds == it & Xinds == ix);

                Tsm(ind) = (Tsb(it) + Tsb(it+1))./2;
                Xsm(ind) = (Xsb(it) + Xsb(it+1))./2;

                % %Calculate mean temperature
                % if it == 1
                %     Tsm(ind) = (T_largest + Tsb(it))./2;
                %     Xsm(ind) = (X_largest + Xsb(it))./2;
                % elseif it == NT
                %     Tsm(ind) = (Tsb(it-1) + Tsb(it))./2;
                %     Xsm(ind) = (Xsb(it-1) + Xsb(it))./2;
                % else
                %     Tsm(ind) = (Tsb(it-1) + Tsb(it))./2;
                %     Xsm(ind) = (Xsb(it-1) + Xsb(it))./2;
                % end

                
                %Tsm(ind) = (Tsb(1:end-1) + Tsb(2:end))./2;
                
    
                


    
            end

        end
        

        
        
    end

end



                    % %Pull all conversions
                    % Xs_cs = zeros(params.Nz, NX);
                    % for ix = 1:NX
                    %     Xs_cs(:, ix) = params.T_z.chars.Xs_cs{im}{it}{ix};
                    % end
                    % Xs_cs_f = Xs_cs(end,:);
                    % 
                    % 
                    % %Calculate average conversion
                    % Xsm(iters) = params.T_z.chars.Xs_cs{im}{it}{ix}(iz) + 1;