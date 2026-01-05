function [T,Y] = TemperatureCharacteristic(t_i, t_f, Ti, y0, params, odeparams, debug)

    %
    im = odeparams.im;
    Xi = odeparams.Xi;
    zi = 0;
    ui = 0;

    %Define initial conditions
    if isempty(y0)
        y0 = [ui; zi; Ti; Xi];
    end

    %Solve ODE system
    t_i = round(t_i, 4);
    options = odeset(RelTol=1e-8,AbsTol=1e-10);
    tspan = [t_i, t_f]; %t_i:0.0001:t_f;
    odefunc = @(t,y) TempOde(t,y,odeparams);
    if any(isnan(tspan)) || any(isnan(y0))
        error('Invalid Input');
    end

    [T,Y] = ode45(odefunc, tspan, y0, options);

     %Debug plots
    if debug & false

        %Plot values versus time
        figure();

        subplot(2,2,1);
        plot(T, Y(:,1), 'r-');
        xlabel('Time (s)'); ylabel('Velocity (m/s)');

        subplot(2,2,2);
        plot(T, Y(:,2), 'k-');
        xlabel('Time (s)'); ylabel('Vertical Pos.');

        subplot(2,2,3);
        plot(T, Y(:,3), 'r-');
        xlabel('Time (s)'); ylabel('Temperature');

        subplot(2,2,4);
        plot(T, Y(:,4), 'k-');
        xlabel('Time (s)'); ylabel('Conversion');


    end

    %Print update
    fprintf('--Characteristic (i_m = %d; T_i = %0.2f K; X_i = %0.3f %%) - Done\n', im, Ti, 100.*Xi);



end

function dydt = TempOde(t,y,params)



    %Pull differential variables
    u = y(1);
    z = y(2);
    T = y(3);
    X = y(4);

    if z < 0
        z = 0;
    end

    %Calculate moles from X
    n_Ar = params.n_Ar_i;
    n_CH4 = params.n_CH4_i * (1-X);
    n_H2 = 2 * params.n_CH4_i * X;
    n_C = 1 * params.n_CH4_i * X;
    n_gas = n_CH4 + n_H2;
    n_tot = n_gas + n_C;

    %Determine local conditions
    T_L = params.T_Lz(z);
    p = params.p_z(z);
    V = (n_gas * params.R * T)/p;
    Db = 2 * ((3 * n_gas * params.R * T)/(4 * pi * p))^(1/3);
    Rb = Db/2;
    As = 4 * pi * (Db/2).^2;
    Ac = As/4; %m2 - cross sectional area
    rho_L = params.rho_L_T(T_L);
    rho_G = (p .* params.methane.mol_mass)./(params.R .* T);
    mu_L = params.liquid.dyn_visc(T);
    sigma = params.liquid.surf_tension(T);
    Eob = ((rho_L - rho_G) .* params.g .* Db.^2)./sigma;

    %Calculate current internal resistance
    n_gas = n_CH4 + n_H2 + n_Ar;
    X_CH4 = n_CH4/n_gas;
    X_H2 = n_H2/n_gas;
    X_Ar = n_Ar/n_gas;
    k = X_CH4 * params.methane.therm_cond(T) + X_H2 * params.hydrogen.therm_cond(T)...
        + X_Ar * params.argon.therm_cond(T);
    Lc = Rb/3;
    Rc = Rb - Lc;
    R_cond = (1/(4 .* pi .* k)) .* (1/Rc - 1/Rb);

    %Calculate current reaction enthalpy
    H_rxn = -6.29E-16 * T^5 + 1.63E-12 * T^4 + 6.82E-9 * T^3 - 3.42E-5*T^2 +0.0469*T + 73.7; 
    H_rxn = 1000*H_rxn; %J/mol-CH4;
    


    %Calculate heat capacity
    C = n_CH4 * params.methane.Cp(T) + n_H2 * params.hydrogen.Cp(T)...
        + n_C * params.carbon.Cp(T) + n_Ar * params.argon.Cp(T);

    %Calculate heat transfer coefficient
    if params.h_conv_const
        h_conv = params.h_conv;
    else
        Nu_conv = 1;
        h_conv = params.h_conv; %UPDATE
    end
    

    %Calculate bubble velocity 
    %ub = CalcVelocity(Db, rho_L, mu_L, sigma, Eob);

    %Calculate reaction rate
    kr = params.k0 .* exp(params.Ea./(params.R .* T));

    %Forces on bubble
    Fb = rho_L * V * params.g;
    Fg = params.m * params.g;
    Cd = CalcCd(u, Db, rho_L, mu_L, sigma, Eob);
    Fd = 0.5 * rho_L * u^2 * Cd * Ac;
    m_vm = (11/16) * rho_L * V;

    %Derivatives
    if z < params.h_reactor

        %Velocity, conversion, and position
        dudt = (Fb - Fd - Fg)/(m_vm + params.m);
        dXdt = kr .* (1 - X);
        dzdt = u;

        %Calculate heat consumption
        dnCH4dt = -dXdt*n_CH4; %mol/s
        Q_rxn = dnCH4dt * H_rxn; %J/s

        %Heat transfer
        if abs(T_L - T) > 1E-4
            if strcmp(params.chars.LC_mode, 'Cond')
                dTdt = (1/(C)) * (T_L - T)/R_cond + Q_rxn/C;
            elseif strcmp(params.chars.LC_mode, 'Conv')
                dTdt = (h_conv * As)/(C) * (T_L - T) + Q_rxn/C;
            end
        else
            dTdt = 0;
        end

    else
        dudt = 0; dTdt = 0; dXdt = 0; dzdt = 0;
    end

    %Catch nans
    if isnan(dudt)
        x = 1;
    end


    %Define output vector
    dydt = [dudt; dzdt; dTdt; dXdt];


end