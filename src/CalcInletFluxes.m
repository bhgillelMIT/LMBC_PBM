%This function calculates the initial fluxes, accounting for 


function N_dot_os = CalcInletFluxes(params)

    if params.heat.active 

        %Handle reaction
        if params.react.active
            error('Not implemented yet. Revisit.')
        end

        N_dot_os = [];

        for im = 1:params.Nms

            %Pull the initial flux 
            init_flux = params.N_dot_o(im);
    
            %Pull initial fraction
            init_fracs = params.T_z.chars.fracs_in{im};
            init_fracs = fliplr(init_fracs);


            %ISSUE HERE - INVERT INDEX FOR CALCULATING FLUXES


            %Log results
            N_dot_os_m = init_fracs .* init_flux;
            N_dot_os = [N_dot_os, N_dot_os_m];
            
    
    
    
    
        end

    else %Simple case of no temp/conv bins
        N_dot_os = params.N_dot_o;

    end


end