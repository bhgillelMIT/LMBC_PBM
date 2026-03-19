function PBM_tests()


%% Create water structure

    water.dyn_visc = @(T) 0.0010016 .* ones(size(T)); %Pa*s - at 20 C
    water.mol_mass = 0.0180153; %kg/mol
    water.surf_tension =  @(T) 0.0728 .*ones(size(T)); 
    water.density = @(T) 998.21 .* ones(size(T)); 
    water.Cp = @(T) 4157 .* ones(size(T)); 
    water.therm_cond = @(T) 0.598 .* ones(size(T)); 
    water.name = 'water';



%% Test 4 - Uniform Binary Breakage

    %Settings
    loadresults = true;
    filepath = 'Data/Solutions/PBM_output_28-Feb-2026_15-38-06.mat';


    %Set inputs
    B0 = 35000;

    inputs = PBM_inputs();
    inputs.src.coalesce_active = false;
    inputs.src.breakage_active = true;
    inputs.src.breakage_model = 'Uniform_Binary'; % 'Hounslow_1988'; 'Scott_1968' 
    inputs.src.breakage_constant_rate = B0;
    inputs.sol.heat.active = false;
    inputs.disc.r_min = 1E-5; %For this case, the unit is m3
    inputs.disc.r_max = 2;
    inputs.disc.Nms = 100;
    inputs.disc.dt = 1;
    inputs.mmesh.type = 'Geometric';
    inputs.mmesh.input = 'Volume'; 
    inputs.sim.t_end = 50;
    inputs.sol.single_layer = true;

    %Run simulation
    if loadresults
        data = load(filepath);
        output = data.output;
        results = PBM_postprocess(output.T, output.Y, output.params);
    else
        results = PBM_v3(inputs);
    end

    %Ode solution
    ys = results.y;
    ts = results.t;
    tsi = ts;
    params = results.params;

    %Will calculate at whatever T is desired 
    Vbs = params.disc.Vbs;
    dVbs = diff(Vbs); dVbs = repmat(dVbs, length(ts), 1);
    yi = ys;
    ys = yi./dVbs;

    %Post process numerical result


    %Calculate numeric density versus time
    Nt = zeros(length(ts), 1);
    for i = 1:length(ts)
        Nt(i) = trapz(yi(i,:));
    end
    Nt = Nt./Nt(1);

    [~,ind6] = min(abs(Nt - 6));
    [~,ind300] = min(abs(Nt - 300));

    %Calculate analytic solution
    vs = logspace(log10(min(Vbs)), log10(0.9999999), 100);
    ns = zeros(3,length(vs));
    ts_comp_inds = 2:100:length(ts)
    ts_comp = ts(ts_comp_inds); %[0; ts(ind6); ts(ind300)];
    for it = 1:length(ts_comp)
        t = ts_comp(it)*B0;
        y_func = @(y) exp(-t.*y.^3)./(y.^2);
        n_func = @(v) exp(-t.^2) .* dirac(v-1) + (6.*t.*v) .* integral(y_func, v, 1); %exp(-t.^2) .* dirac(v-1) + 6.*t.*v .* integral(y_func, v, 1);
        %Iterate through each v
        for iv = 1:length(vs)
            v = vs(iv);
            ns(it, iv) = (1/7.7545)*n_func(v);
        end

    end


    %Calculate number
    Nt_analytical = zeros(1, length(ts_comp));
    for i = 1:length(ts_comp)
        Nt_analytical(i) = trapz(ns(i,:));
    end
    Nt_analytical = Nt_analytical./Nt_analytical(1);

    figure();
    plot(ts_comp, Nt_analytical, 'k-'); hold on;
    plot(ts_comp, Nt(ts_comp_inds), 'ko');
    
    %ns(1,:) = n_func(vs, 0);
    %ns(2,:) = n_func(vs, ts(ind6));
    %ns(3,:) = n_func(vs, ts(ind300));


    %Plot solution at several times
    figure();



    
    loglog(params.Vms, ys(1,:), 'k.'); hold on;
    loglog(params.Vms, ys(ind6,:), 'ko');
    loglog(params.Vms, ys(ind300,:), 'k+');

    loglog(vs, ns(2,:));
    loglog(vs, ns(3,:));
    %loglog(params.Vms, n_constant(end,:), 'k-');
    %loglog(params.Vms, n_constant(2,:), 'k--');
    %loglog(params.Vms, n_constant(1,:), 'k:')
    %loglog(params.Vms, Fvals, 'k--');
    ylim([1E-5, 10000]); xlim([inputs.disc.r_min, inputs.disc.r_max])
    xlabel('Particle Volume (m^3)'); ylabel('Number density (1/m^3)');
    axis square; grid on; 

    % subplot(1,2,2);
    % semilogx(xs, n_constant_normal);
    



%% Test 3 - Constant Coalescence Numerical Test 


    %Specify output folder
    output_folder = 'Data/Studies/Analytical/';


    %Set inputs
    C = 1000;
    N0 = 1;
    v0 = 0.01;
    inputs = PBM_inputs();
    inputs.src.coalesce_active = true;
    inputs.src.breakage_active = false;
    inputs.src.coalesce_model = 'Scott_1968'; % 'Hounslow_1988'; 'Scott_1968' 
    inputs.src.coalesce_constant_rate = C;
    inputs.sol.heat.active = false;
    inputs.disc.r_min = 1E-5; %For this case, the unit is m3
    inputs.disc.r_max = 10000;
    inputs.disc.Nms = 150;
    inputs.disc.dt = 0.0025;
    inputs.mmesh.type = 'Geometric';
    inputs.mmesh.input = 'Volume'; 
    inputs.sim.t_end = 5;
    inputs.sol.single_layer = true;



    %Run simulation
    results = PBM_v3(inputs);


    %Ode solution
    ys = results.y;
    ts = results.t;
    params = results.params;

    %Will calculate at whatever T is desired 
    Vbs = params.disc.Vbs;
    dVbs = diff(Vbs); dVbs = repmat(dVbs, length(ts), 1);
    yi = ys;
    ys = yi./dVbs;
    %ys = yi;

    %Calculate numeric density versus time
    Nt = zeros(length(ts), 1);
    for i = 1:length(ts)
        Nt(i) = trapz(params.Vms, ys(i,:));
    end
    Nt = Nt./Nt(1);



    %Calculate analytical soluti on
    ts_ref = [0.00000000001, ts(6), ts(484)];






    switch inputs.src.coalesce_model
        case 'Hounslow_1988'
            [n_constant, phi_constant] =  Hounslow_1988(ts_ref, N0, v0, C, params); % Scott_1968(ts_ref, N0, 2*v0, C, params); %Hounslow_1988(ts_ref, N0, v0, C, params); %

        case 'Scott_1968'
            [n_constant, phi_constant] =  Scott_1968(ts_ref, N0, v0, C, params); %Hounslow_1988(ts_ref, N0, v0, C, params); %
    end

    %Calculate normalized curve
    vs = params.Vms;
    xs = vs/v0;
    n_constant_normal = xs.^2 .* phi_constant;

    %Convert to normal density



    

    %Calculate number over time
    Nt = zeros(length(ts), 1);
    for it = 1:length(Nt)
        


    end



    %Plot solution at several times
    figure();

    F_func = @(v) (N0/v0) .* exp(-v./v0);
    Fvals = F_func(params.Vms);

    
    loglog(params.Vms, ys(1,:), 'k.'); hold on;
    loglog(params.Vms, ys(6,:), 'ko');
    loglog(params.Vms, ys(484,:), 'k+');
    loglog(params.Vms, n_constant(end,:), 'k-');
    loglog(params.Vms, n_constant(2,:), 'k--');
    loglog(params.Vms, n_constant(1,:), 'k:')
    %loglog(params.Vms, Fvals, 'k--');
    ylim([1E-15, 1000]); xlim([inputs.disc.r_min, inputs.disc.r_max])
    xlabel('Particle Volume (m^3)'); ylabel('Number density (1/m^3)');
    axis square; grid on; 

    % subplot(1,2,2);
    % semilogx(xs, n_constant_normal);
    


   %N_t = @(t) N0 * exp()
    %tau = b * N0 * v0 * t;




    %Save results

%% Test 1 - Comparison to Wang et al 1D

    %Specify output folder
    output_folder = 'Data/Studies/Flow Rate/';

    %Test conditions
    u_spfs = [0.01, 0.03, 0.05, 0.12, 0.16];
    alpha_gs = [0.038 0.109, 0.154, 0.198, 0.225];
    epss = 10 .* u_spfs;

    %Create storage array
    results = cell(1, length(u_spfs));
    ys_out = cell(1, length(u_spfs));



    %Iterate through cases
    for ic = 1:length(u_spfs)

        %Modify inputs
        inputs = PBM_inputs();
        inputs.reactor.liquid = water;
        inputs.reactor.u_spf_orifice = u_spfs(ic); %m/s
        inputs.reactor.T_gas_i_mu = 20 + 273.15;
        inputs.reactor.T = 20 + 273.15;
        inputs.reactor.T_min_i = 0 + 273.15;
        inputs.src.eps_manual = epss(ic);
        inputs.src.solve_eps = false;
        inputs.src.alphag_manual = alpha_gs(ic);
        inputs.src.solve_alphag = false;
        inputs.sol.heat.active = false;
        inputs.sol.react.active = false;
        inputs.sol.solve_ub = false;
        inputs.sol.ub_manual = inputs.reactor.u_spf_orifice/inputs.src.alphag_manual; 
        inputs.reactor.H = 0.3;
        inputs.sim.t_end = 3;
        inputs.sol.src_delay = 0;
        inputs.sol.break_file = 'break_water_Nd-20_Nz-1_Ne-20_TL-293.mat';

        %Run PBM
        results{ic} = PBM_v3(inputs);
        ys_out{ic} = results{ic}.y;
        ys_final{ic} = results{ic}.y(end,:);
        
    end

    %Save data
    outstruct.results = results; outstruct.ys_out = ys_out;
    outstruct.ys_final = ys_final;
    outfilename = sprintf('outstruct.mat');
    save(outfilename, 'outstruct');

    
%% Plot results

    %Load data
    data = load('outstruct.mat'); data = data.outstruct;
    params = load('Data/Studies/Flow Rate/flow_rate_study_params.mat'); 
    params = params.params;

    figure();

    for i = 1:5

        %Calculate equivalent diameters
        p = params.p_surf;
        T_bar = params.T_liq;
        ns_gas = params.nms;
        Vs = (ns_gas * params.R * T_bar)/p;
        ds = (6.*Vs/pi).^(1/3);


        y_final = data.ys_final{i}(1:30);

        
        subplot(1,2,1);
        semilogx(1000.*params.mms, y_final./max(y_final), 'LineWidth', 2); hold on;

        subplot(1,2,2);
        plot(ds, y_final./max(y_final), 'LineWidth', 2); hold  on;

    end

    subplot(1,2,1);
    xlabel('Bubble Mass (g)'); ylabel('Normalized BSD');
    grid on; grid minor; axis square;
    legend()
    set(gca, 'FontSize', 18, 'FontWeight', 'bold');

    subplot(1,2,2);
    xlabel('Bubble Diameter (m)'); ylabel('Normalized BSD');
    grid on; axis square;
    set(gca, 'FontSize', 18, 'FontWeight', 'bold');

    

    x=1

%% Test 2 - Flow rate study

    u_spfs = 0.01:0.01:0.05;

    inputs = PBM_inputs();

    for iu = 1:length(u_spfs)

        %Modify inputs
        inputs = PBM_inputs();
        inputs.reactor.liquid = water;
        inputs.reactor.u_spf_orifice = u_spfs(iu); %m/s
        inputs.reactor.T_gas_i_mu = 20 + 273.15;
        inputs.reactor.T = 20 + 273.15;
        inputs.reactor.T_min_i = 0 + 273.15;
        inputs.sol.heat.active = false;
        inputs.sol.react.active = false;
        inputs.reactor.H = 1;
        inputs.sim.t_end = 5;
        inputs.sol.src_delay = 3;

        %Run PBM
        results = PBM_v3(inputs);


    end








end



function [n_constant, phi_constant_out] = Scott_1968(ts, N0, v0, C, params)

    


    %Inouts
    nu = 1;
    
    N0 = N0; % number of droplets initially


    %Calculate the true mean
    n_func = @(v) (N0/v0) .* (v/v0) .* exp(-v/v0);
    v_func = @(v) v .* n_func(v);
    v_bar = integral(v_func, 1E-9, 1E5)/integral(n_func, 1E-9, 1E5);
    v0 = v_bar;

    N_sum_vals = 30;

    %Calculate input values 
    vs = params.Vms;
    xs = vs/v0;
    Ts = C.*N0.*ts;

    %Create output storage vector
    n_constant = zeros(length(ts), length(params.Vms));
    phi_constant_out = n_constant;
    Nt_analytic = zeros(length(ts));

    for i = 1:length(ts)

        %Calculate current non-dim time
        t = ts(i);
        T = C.*N0.*t;

        %Calculate sum
        Nt_analytic(i) = 2.*N0./(T+2);
        phi_sum_func = @(x,T,k) (((x.*(nu+1)).^((nu+1).*(k+1)))./(gamma((nu+1).*(k+1)))) .* (T./(T+2)).^k;


        % %Calculate sum at this time
        % phi_sum_vals = zeros(N_sum_vals, length(xs));
        % for k = 0:(N_sum_vals-1)
        %     phi_sum_vals(k+1,:) = phi_sum_func(xs, T, k);
        % end
        % phi_sum_vals(isinf(phi_sum_vals)) = 0;
        % phi_sum_vals(isnan(phi_sum_vals)) = 0;
        % 
        % phi_sum = sum(phi_sum_vals);

        %Calculate non-dim spectrum

        sinh_ins = 2.*xs.*sqrt(T./(T+2));
        exp_ins = -2.*xs;
        phi_constant = (8 .* exp(exp_ins) .* sinh(sinh_ins))./...
            (sqrt(T) .* (T+2).^1.5);
        
        % phi_constant = zeros(1,length(params.Vms));
        % for ix = 1:length(params.Vms)
        %     x = xs(ix);
        %     phi_sum_vals = zeros(1, N_sum_vals);
        %     for k = 0:(N_sum_vals-1)
        %         phi_sum_vals(k+1) = phi_sum_func(x, T, k);
        %     end
        %     phi_sum_vals(isinf(phi_sum_vals)) = 0;
        %     phi_sum_vals(isnan(phi_sum_vals)) = 0;
        %     phi_sum = sum(phi_sum_vals);
        %     phi_constant(ix) = (4.*exp(-(nu+1).*x))./(x.*(T+2).^2) .* phi_sum; %4 .* exp(-2.*xs./(T+2))./((T+2).^2); %%(8 .* exp(-2.*xs) .* sinh(2.*xs.*sqrt(T./(T+2))))./(sqrt(T) .* (T+2).^(1.5));  %(4.*exp(-(v+1).*xs))./(xs.*(T+2).^2) .* sum(phi_sum_vals); %4 .* exp(-2.*xs./(T+2))./((T+2).^2); %
        % end
        phi_constant(isnan(phi_constant)) = 0;
        phi_constant_out(i,:) = phi_constant;

        %Convert to actual spectrum
        
        n_constant(i, :) = (N0/v0) .* phi_constant;
        

    end
    


end



function [n_constant, phi_constant_out] = Hounslow_1988(ts, N0, v0, C, params)
 

    %Create output storage vector
    v0 = v0*2;
    N0 = 0.2*N0;
    n_constant = zeros(length(ts), length(params.Vms));
    n_tilde = n_constant;
    Nt_analytic = zeros(length(ts));

    vs = params.Vms;

    for i = 1:length(ts)

        t = ts(i);
        tau = N0*C*t;
        L_tilde = (vs./v0).^(1/3);
        n_tilde(i,:) = (2/(2+tau)) .* (exp(-(2*L_tilde.^3)/(tau + 2)) - exp(-4.*L_tilde.^3./(tau + 2)));
    

        n_prime = (4.*N0)./(v0 .* (tau + 2).^2) .* exp(-(2 .* vs./v0)./(tau + 2));

        n_constant(i,:) = n_tilde(i,:) * N0; %n_prime; %_tilde(i,:) * N0; %n_tilde(i,:) * N0;
    end

    phi_constant_out = 1;



end
