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
        merged = params.T_z.chars.merged{im};
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