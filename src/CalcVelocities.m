function ubs = CalcVelocities(Ds, Ts, liquid, p_orifice, fsolve_opts)

    %Physical constants
    g = 9.81;
    R = 8.3145;
    M_CH4 = 0.016;

    %Iterate through all combinations
    for it = 1:length(Ts)

        %Pull temperature
        T = Ts(it);

        %Calculate properties
        rho_G = (p_orifice .* M_CH4)./(R .* T);
        rho_L = liquid.density(T);
        mu_L = liquid.dyn_visc(T);
        sigma = liquid.surf_tension(T);


        %Iterate through diameters
        for id = 1:length(Ds)

            %
            Db = Ds(id);

            %Calculate eotvos number
            Eob = ((rho_L - rho_G) .* g .* Db.^2)./sigma;
            
            %Calculate velocity
            ub = CalcVelocity(Db, rho_L, mu_L, sigma, Eob, fsolve_opts);

            %Log result
            ubs(id, it) = ub;

        end
    end


end