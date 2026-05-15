%This function calculates the initial fluxes, accounting for 


function N_dot_os = CalcInletFluxes(params)

    if params.heat.active & params.react.active

        N_dot_os = [];

        %Pull conversions of the largest case
        Tsb_largest = interp1(params.T_z.chars.largest.zs, params.T_z.chars.largest.Ts, params.mesh.yy(:,1));
        Xsb_largest = interp1(params.T_z.chars.largest.zs, params.T_z.chars.largest.Xs, params.mesh.yy(:,1));
        

        %Iterate through each representative size
        for im = 1:params.Nms

            %Pull the initial flux 
            init_flux = params.N_dot_o(im);
    
            %Pull initial fraction for temperatures
            init_fracs = params.T_z.chars.fracs_in{im};
            init_fracs = fliplr(init_fracs);

            %Calculate flux into these temperature bins
            N_dot_os_m_T = init_fracs .* init_flux;

            %Pull trailing inds
            itrails = params.T_z.chars.itrails{im};

            %Allocate all the mass to the first conversion bin
            N_dot_os_m = [];
            itrail = params.T_z.chars.itrails{im};
            for it = 1:length(init_fracs)

                if any(it == itrails) %Should be assigned to the largest bin (nothing will be assigned though)
                    init_ind = params.chars.N_Xis;
                    N_dot_os_m_T_X(init_ind) = N_dot_os_m_T(it);
                else %Should be assigned to the bin that isn't a trailing conv char, but starts at Xi = 0
    
                    %Identify the conversion characteristic that is the
                    N_dot_os_m_T_X = zeros(1, params.chars.N_Xis);
                    
                    
                    Xsb = zeros(params.Nz+1, params.chars.N_Xis);
                    Xsc = zeros(params.Nz, params.chars.N_Xis);
                    trailing_conv = false(1, params.chars.N_Xis);
                    for ix = 1:params.chars.N_Xis
                        
                        %Pull the conversions
                        Xsb(:,ix) = params.T_z.chars.Xs_bs{im}{it}{ix};
                        Xsc(:,ix) = params.T_z.chars.Xs_cs{im}{it}{ix};

                        %Determine if this is a trailing conversion bin
                        matches_largest = Xsb(:,ix) == Xsb_largest;       
                        trailing_conv(ix) = any(matches_largest(2:end));

    
                    end

                    %Identify trailing conv characteristics
                    nontrailing_conv = ~trailing_conv;
                    correct_Xi = Xsb(1,:) == params.reactor.X_i;
                    valid_ind = find(nontrailing_conv & correct_Xi);
                    if length(valid_ind) == 1
                        N_dot_os_m_T_X(valid_ind) = N_dot_os_m_T(it);
                    else
                        error('Too few or too many valid starting reaction bins. Programming error.');
                    end


                    %N_dot_os_m_T_X(1) = N_dot_os_m_T(it);

                end

                %Set to zero if it is a trailing bin 
                if any(it == itrail)
                    N_dot_os_m_T_X = zeros(size(N_dot_os_m_T_X));
                end


                %Add to group
                N_dot_os_m = [N_dot_os_m, N_dot_os_m_T_X];
            end


            


            %Log results
            
            N_dot_os = [N_dot_os, N_dot_os_m];


        end


    elseif params.heat.active 

        
        

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