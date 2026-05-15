function params = CalcLocalProperties(y, params, local)

    %Handle inputs
    if nargin < 3
        local = false;
    end

    %Make params global if not local
    % if ~local
    %     global params
    % end

    %Resolve local gas velocity
    %params.ug = 0.3; %TEMPORARY

    %Determien iz inds to use
    if params.sol.sep_layer & strcmp(params.sol.type, 'segregated')
        iz_inds = params.iz;
    else
        iz_inds = 1:params.Nz;
    end

    %Iterate through spatial cells
    for iz = iz_inds

        %Pull current z
        z = params.zms(iz);

        cellinds = ((iz-1).*params.Nms + 1):1:(iz.*params.Nms);

        %Pull numeric densities in this spatial cell
        if all(y == y(1)) %Source test case
            Ns_cell = y(cellinds);
        else
            Ns_cell = params.Ns_m(cellinds); %y(cellinds);
        end


        %Calculate gas holdup
        V_cell = params.Vcells(iz);
        V_gas = 0;
        V_dot = zeros(1, params.Nms);
        for ix = 1:params.Nms

            %Calculate current mean diameters
            params.n_mu(iz, ix) = params.nms(ix) .* params.X_Ar + params.nms(ix) .* (1 + params.X_mu(iz,ix));
            params.V_mu(iz, ix) = (params.n_mu(iz, ix) .* params.R .* params.T_mu(iz,ix))./params.p_func(z);
            params.d_mu(iz, ix) = ((6 .* params.V_mu(iz, ix))./pi).^(1/3);

            if any(~isreal(params.d_mu))
                
            end


            %Calculate current bubble velocities
            if params.sol.solve_ub
                if ~isreal(params.d_mu(iz,ix))
                    params.uz_mu(iz, ix) = 0.3; %params.ubs.funcs{iz}();
                else
                    params.uz_mu(iz, ix) = params.ubs.funcs{iz}(params.d_mu(iz, ix));
                end

                %params.uz_mu(iz, ix) = params.ubs.funcs{iz}(params.d_mu(iz, ix));
            else
                if length(params.sol.ub_manual) == params.Nms
                    params.uz_mu(iz, ix) = params.sol.ub_manual(ix);
                else
                    params.uz_mu(iz, ix) = params.sol.ub_manual(1);
                end
            end

            %Calculate gas flow rates
            V_dot(ix) = Ns_cell(ix) * V_cell * params.V_mu(iz, ix) * (params.uz_mu(iz,ix)./params.mesh.hs(iz));
            V_gas = V_gas + Ns_cell(ix) * V_cell * params.V_mu(iz, ix);

        end

        %Calculate terms at boundaries
        for ix = 1:(params.Nms + 1)

            %Calculate current mean diameter
            if ix < (params.Nms + 1) 
                params.nb_mu(iz, ix) = params.nbs(ix) .* params.X_Ar + params.nbs(ix) .* (1 + params.X_mu(iz,ix));
                params.Vb_mu(iz, ix) = (params.nb_mu(iz, ix) .* params.R .* params.T_mu(iz,ix))./params.p_func(z);
            else
                params.nb_mu(iz, ix) = params.nbs(ix) .* params.X_Ar + params.nbs(ix) .* (1 + params.X_mu(iz,ix-1));
                params.Vb_mu(iz, ix) = (params.nb_mu(iz, ix) .* params.R .* params.T_mu(iz,ix-1))./params.p_func(z);
            end
            params.db_mu(iz, ix) = ((6 .* params.Vb_mu(iz, ix))./pi).^(1/3);

            %Calculate differences
            if ix > 1
                params.dbds(iz, ix-1) = params.db_mu(iz, ix) - params.db_mu(iz, ix-1);
            end

            
     

        end



        V_dot_total = sum(V_dot); %m^3/s

        %Calculate gas holdup
        if V_gas < 0
            V_gas = 0;
        end
        params.alpha_g(iz) = V_gas/V_cell;

        %Calculate gas superficial velocity in each cell
        params.ug(iz) = V_dot_total/params.reactor.Ac;


        %Calculate turbulent energy dissipation
        if params.src.solve_eps
            params.turb.eps(iz) = params.g * params.ug(iz)+1E-12;
        else
            params.turb.eps(iz) = params.src.eps_manual;
        end
        params.turb.k(iz) = 0.1;
            %turb.eps_G(iz) = 1;

        %Impose manual gas holdup if specifified 
        if ~params.src.solve_alphag
            params.alpha_g(iz) = params.src.alphag_manual;
        end
         
        %Check if value is reasonable
        if params.alpha_g(iz) > params.alpha_g_max
            params.alpha_g(iz) = params.alpha_g_max - 1E-6;
            x = 1;
        end
        
    end

    %Check outputs for validity
    if any(isnan(params.turb.eps) | isnan(params.alpha_g))
        warning('NaNs in local properites.')
    elseif params.turb.eps < 0
        warning('Negative turbulent eps.')
    elseif any(params.alpha_g < 0)
        warning('Negative gas holdup.');
    end

end