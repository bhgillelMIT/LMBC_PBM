function PBM_flowratesweep()


%% Create water structure

    water.dyn_visc = @(T) 0.0010016 .* ones(size(T)); %Pa*s - at 20 C
    water.mol_mass = 0.0180153; %kg/mol
    water.surf_tension =  @(T) 0.0728 .*ones(size(T)); 
    water.density = @(T) 998.21 .* ones(size(T)); 
    water.Cp = @(T) 4157 .* ones(size(T)); 
    water.therm_cond = @(T) 0.598 .* ones(size(T)); 
    water.name = 'water';

%% Test 1 - Comparison to Wang et al 1D

    %Specify output folder
    output_folder = 'Data/Studies/Flow Rate/';

    %Test conditions
    u_spfs = [0.01, 0.02, 0.03, 0.04, 0.05]; %[0.01, 0.03, 0.05, 0.12, 0.16];
    alpha_gs = [0.038 0.109, 0.154, 0.198, 0.225];
    epss = 10 .* u_spfs;

    %Create storage array
    results = cell(1, length(u_spfs));
    ys_out = cell(1, length(u_spfs));



    %Iterate through cases
    for ic = 1:length(u_spfs)

        %Modify inputs
        inputs = PBM_inputs();
        inputs.disc.Nms = 60;
        inputs.disc.mesh_hybrid_cells = round(inputs.disc.Nms .* 0.5);  
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
        inputs.sol.single_layer = true;
        inputs.sol.ub_manual = inputs.reactor.u_spf_orifice/inputs.src.alphag_manual; 
        inputs.reactor.H = 0.3;
        inputs.src.breakage_active = true;
        inputs.sim.t_end =1.5;
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
    timestamp = datestr(now, 'yyyy-mm-dd-HH-MM-SS');
    outfilename = sprintf('flowrate_%s.mat', timestamp);
    outfilepath = fullfile(output_folder, outfilename);
    save(outfilepath, 'outstruct');

    
%% Plot results

    %Clean up
    close all

    %Specify plot colors
    colors = {[0.85, 0.3250, 0.0980];
              [0.929, 0.694, 0.125];
             [0.494, 0.184, 0.556];
              [0.466, 0.674, 0.188];
               [0.301, 0.745, 0.933]};

    %Load data
    data = load('Data/Studies/Flow Rate/flowrate_2026-03-19-11-55-22.mat'); 
    data = data.outstruct;
    %params = load('Data/Studies/Flow Rate/flow_rate_study_params.mat'); 
    %params = params.params;
    params = data.results{1}.params;

    initfig = figure();
    normfig = figure();
    actfig = figure();
    pdffig = figure();

    for i = 1:5

        %Calculate equivalent diameters
        p = params.p_surf;
        T_bar = params.T_liq;
        ns_gas = params.nms;
        Vs = (ns_gas * params.R * T_bar)/p;
        ds = (6.*Vs/pi).^(1/3);

        ds = ds(2:end);


        mms = params.mms(2:end);

        %Pull and plot initial distribution
        y_initial = data.ys_out{i}(1, 2:end);
        figure(initfig);
        subplot(1,2,1);
        semilogx(1000.*mms, y_initial, 'LineWidth', 2); hold on;
        subplot(1,2,2);
        plot(ds, y_initial, 'LineWidth', 2); hold  on;





        %Plot final result
        y_final = data.ys_final{i}(2:end);

        figure(actfig);
        subplot(1,2,1);
        semilogx(1000.*mms, y_initial, 'k--', 'LineWidth', 2, 'Color', colors{i}); hold on;
        semilogx(1000.*mms, y_final, 'k-', 'LineWidth', 2, 'Color', colors{i}); hold on;
        subplot(1,2,2);
        plot(ds, y_initial, 'k--', 'LineWidth', 2, 'Color', colors{i}); hold on;
        plot(ds, y_final, 'k-', 'LineWidth', 2, 'Color', colors{i}); hold  on;

        figure(normfig);
        subplot(1,2,1);
        semilogx(1000.*mms, y_final./max(y_final), 'LineWidth', 2); hold on;
        subplot(1,2,2);
        plot(ds, y_final./max(y_final), 'LineWidth', 2); hold  on;

        %Plot probability density function in terms of volume
        figure(pdffig);
        Vbs = Vs(2:end) .* y_final;
        V_total = sum(Vbs);

        plot(ds, Vbs./V_total, 'LineWidth', 2, 'Color', colors{i}); hold on;

    end

    %Create legend for actual figure
    figure(actfig);
    subplot(1,2,2);
    h = zeros(1, 7);
    h(1) = plot(nan, nan, 'k--', 'LineWidth', 2); hold on;
    h(2) = plot(nan, nan, 'k-', 'LineWidth', 2); hold on;
    h(3) = plot(nan, nan, 'k-', 'LineWidth', 2, 'Color', colors{1}); hold on;
    h(4) = plot(nan, nan, 'k-', 'LineWidth', 2, 'Color', colors{2}); hold on;
    h(5) = plot(nan, nan, 'k-', 'LineWidth', 2, 'Color', colors{3}); hold on;
    h(6) = plot(nan, nan, 'k-', 'LineWidth', 2, 'Color', colors{4}); hold on;
    h(7) = plot(nan, nan, 'k-', 'LineWidth', 2, 'Color', colors{5}); hold on;

    %Format actual figure
    subplot(1,2,1);
    %legend(h, 'initial', 'Final', 'u = 0.01 m/s', 'u = 0.03 m/s', 'u = 0.05 m/s', 'u = 0.12 m/s', 'u = 0.16 m/s' ,'location', 'northwest');
    xlabel('Bubble Mass (g)'); ylabel('Actual BSD');
    grid on; grid minor; axis square;
    set(gca, 'FontSize', 18, 'FontWeight', 'bold');

    subplot(1,2,2);
    xlabel('Bubble Diameter (m)'); ylabel('Actual BSD');
    grid on; axis square;
    set(gca, 'FontSize', 18, 'FontWeight', 'bold');
    legend(h, 'initial', 'Final', 'u = 0.01 m/s', 'u = 0.02 m/s', 'u = 0.03 m/s', 'u = 0.04 m/s', 'u = 0.05 m/s' ,'location', 'northwest');


    %Format normalized figure
    figure(normfig);
    subplot(1,2,1);
    xlabel('Bubble Mass (g)'); ylabel('Normalized BSD');
    grid on; grid minor; axis square;
    legend('initial','u = 0.01 m/s', 'u = 0.03 m/s', 'u = 0.05 m/s', 'u = 0.12 m/s', 'u = 0.16 m/s');
    set(gca, 'FontSize', 18, 'FontWeight', 'bold');

    subplot(1,2,2);
    xlabel('Bubble Diameter (m)'); ylabel('Normalized BSD');
    grid on; axis square;
    
    set(gca, 'FontSize', 18, 'FontWeight', 'bold');

    %Format PDF figure
    figure(pdffig);
    xlabel('Bubble Diameter (m)'); ylabel('PDF (volume) of d_b 1/m');
    grid on; axis square;
    legend('u = 0.01 m/s', 'u = 0.03, m/s', 'u = 0.05 m/s', 'u = 0.12 m/s', 'u = 0.16 m/s');
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