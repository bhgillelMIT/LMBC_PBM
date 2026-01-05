function PBM_tests()


%% Create water structure

    water.dyn_visc = @(T) 0.0010016 .* ones(size(T)); %Pa*s - at 20 C
    water.mol_mass = 0.0180153; %kg/mol
    water.surf_tension =  @(T) 0.0728 .*ones(size(T)); 
    water.density = @(T) 998.21 .* ones(size(T)); 
    water.Cp = @(T) 4157 .* ones(size(T)); 
    water.therm_cond = @(T) 0.598 .* ones(size(T)); 
    water.name = 'water';



%% Test 1 - Comparison to Wang et al 1D

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
    save('outstruct.mat', 'outstruct');

    
%% Plot results

    %Load data
    data = load('outstruct.mat'); data = data.outstruct;
    params = load('flow_rate_study_params.mat'); 
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
    grid on; grid minor; axis square;
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