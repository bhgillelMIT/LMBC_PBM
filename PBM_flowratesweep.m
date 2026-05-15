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

if true
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
    parfor ic = 1:length(u_spfs)

        %Initialize inputs and specify coalescence/breakage kernel
        inputs = PBM_inputs();
        inputs.src.coalesce_model = 'Wang_2005';
        inputs.src.breakage_model = 'Wang_2005';

        %Modify inputs
        
        inputs.disc.Nms = 60;
        inputs.disc.mesh_hybrid_cells = round(inputs.disc.Nms .* 0.5);  
        inputs.reactor.liquid = water;
        inputs.reactor.u_spf_orifice = u_spfs(ic); %m/s
        inputs.reactor.T_gas_i_mu = 20 + 273.15;
        inputs.reactor.T = 20 + 273.15;
        inputs.reactor.T_orifice = inputs.reactor.T;
        inputs.mmesh.type = 'Linear';
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
        inputs.sim.t_end = 3;
        inputs.sol.src_delay = 0;
        inputs.sol.break_file = 'break_water_Nd-20_Nz-1_Ne-20_TL-293.mat';
        inputs.sol.postprocess = false; %Disabled to prevent figures 
        inputs.src.c_wake = 0.3;

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

end
    
%% Plot results

if true

    %Clean up
    close all

    %Specify plot colors
    colors = {[0.85, 0.3250, 0.0980];
              [0.929, 0.694, 0.125];
             [0.494, 0.184, 0.556];
              [0.466, 0.674, 0.188];
               [0.301, 0.745, 0.933]};

    %Load data
    data = load('Data/Studies/Flow Rate/flowrate_2026-05-10-19-02-01.mat'); 
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
        V_PDF = Vbs./V_total;
        V_PDF_sum = sum(V_PDF);
        V_PDF_trapz = trapz(1000*ds, V_PDF);
        V_PDF = V_PDF/V_PDF_trapz;

        subplot(2,1,1);
        plot(ds, y_final, 'LineWidth', 2, 'Color', colors{i}); hold on;

        subplot(2,1,2);
        plot(ds, V_PDF, 'LineWidth', 2, 'Color', colors{i}); hold on;

    
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

    subplot(2,1,1);
    xlabel('Bubble Diameter (m)'); ylabel('Numeric Density (1/m3)');
    grid on; axis square;
    legend('u = 0.01 m/s', 'u = 0.03, m/s', 'u = 0.05 m/s', 'u = 0.12 m/s', 'u = 0.16 m/s');
    set(gca, 'FontSize', 14, 'FontWeight', 'bold');

    subplot(2,1,2);
    xlabel('Bubble Diameter (m)'); ylabel('PDF (volume) of d_b (1/m)');
    grid on; axis square;
    legend('u = 0.01 m/s', 'u = 0.03, m/s', 'u = 0.05 m/s', 'u = 0.12 m/s', 'u = 0.16 m/s');
    set(gca, 'FontSize', 14, 'FontWeight', 'bold');
    

end

%% Test 2 - Full Flow rate study - Understanding evolution of BSD, and gas holdup - ISOTHERMAL NON-REACTING CASE

    u_spfs = [0.0025:0.0025:0.01, 0.02:0.01:0.15];

   
    results_all = cell(1, length(u_spfs));



    %Iterate through cases 
    parfor iu = 1:length(u_spfs)

        %Load Inputs
        inputs = PBM_inputs();

        %Modify Inputs
        inputs.disc.Nms = 100;
        inputs.disc.Nz = 8;
        inputs.disc.r_min = 0.25 .* 0.001; %m - smallest size of bubble to consider
        inputs.disc.r_max = 25 .* 0.001;
        inputs.inlet.d.mu_i = 0.005;
        inputs.inlet.d.std_i = 0.00075;
        inputs.reactor.liquid = water;
        inputs.reactor.u_spf_orifice = u_spfs(iu); %m/s
        inputs.reactor.T_gas_i_mu = 20 + 273.15;
        inputs.reactor.T = 20 + 273.15;
        inputs.reactor.T_min_i = 0 + 273.15;
        inputs.sol.heat.active = false; 
        inputs.sol.react.active = false;
        inputs.sim.t_end = 3;
        inputs.sol.src_delay = 0;
        inputs.sol.break_file = 'break_water_Nd-20_Nz-1_Ne-20_TL-293.mat';
        inputs.sol.postprocess = false; %Disabled to prevent figures 
        inputs.src.c_wake = 1.2;
        inputs.reactor.H = 0.3;
        inputs.src.coalesce_model = 'Wang_2005';
        inputs.src.breakage_model = 'Wang_2005';
        

        %Run PBM
        results = PBM_v3(inputs);

        %Compile results
        results_all{iu} = results;

    end




    %Save the file 
    filepath = 'Data/Studies/Flow Rate/';
    filename = sprintf('flow_rate_sweep_full_%s.mat', datestr(now, 'yyyy-mm-dd-HH-MM-SS'));
    save(fullfile(filepath, filename), 'results_all');

%% Test 2 - plots

    %Plot settings 
    lw = 2;
    fs = 16;

    %Load file
    results_all = load('Data/Studies/Flow Rate/flow_rate_sweep_full_2026-05-13-16-54-28.mat');
    results_all = results_all.results_all;

    %Create storage arrays
    u_spfs = zeros(1, length(results_all));
    alphags = zeros(1, length(results_all));
    d_bars = zeros(1, length(results_all));


    %Tbars_final = zeros(1, length(results_all));
    %Xbars_final = zeros(1, length(results_all));
    BSDs = cell(1, length(results_all));
    BSDs_norm = cell(1, length(results_all));
    dbsf = cell(1, length(results_all));
    V_PDFs = cell(1, length(results_all));

    %Create figure for BSD evolution plot
    BSDfig = figure();
    BSDcmap = nebula(length(results_all));

    %Iterate through cases 
    for ir = 1:length(results_all)

        %Pull results and params
        results = results_all{ir};
        params = results.params;

        %Run through post processing
        results_processed = PBM_postprocess(results.t, results.y, results.params);

        %Extract and plot BSD evolution, gas holdup evolution, final temperature profile, and final conversion profile for each case
        u_spfs(ir) = params.u_spf_orifice;
        alphags(ir) = results_processed.alphag_avg;
        d_bars(ir) = results_processed.dsauter(end);
        BSDs{ir} = results_processed.mdists(end,:);
        BSDs_norm{ir} = BSDs{ir}./trapz(results_processed.dbsout(end,:), BSDs{ir});
        dbsf{ir} = results_processed.dbsout(end,:);

        %Calculate mass/volume probability distribution
        Vbs = params.Vms(1:end) .* BSDs{ir};
        V_total = sum(Vbs);
        V_PDF = Vbs./V_total;
        V_PDF_sum = sum(V_PDF);
        V_PDF_trapz = trapz(1000*dbsf{ir}, V_PDF);
        V_PDF = V_PDF/V_PDF_trapz;
        V_PDFs{ir} = V_PDF;

        %Determien color for 

        %Plot the BSD evolution for this case
        figure(BSDfig);
        plot(dbsf{ir}, BSDs{ir}, 'b-','LineWidth', lw, 'Color', BSDcmap(ir, :)); hold on; hold on;



    end

    %Format BSD figure
    figure(BSDfig);
    xlabel('Bubble Diameter (m)', 'FontSize', fs);
    ylabel('Bubble Size Distribution (1/m3)', 'FontSize', fs);
    grid on; grid minor; axis square;


    %Plot Evolution of bubble diameter and gas holdup
    figure();

    subplot(1,2,1);
    plot(100*u_spfs, 100*alphags, 'k-', 'LineWidth' ,lw);
    xlabel('Superficial Gas Velocity (cm/s)', 'FontSize', fs, 'FontWeight', 'bold');
    ylabel('Gas Holdup (%)', 'FontSize', fs, 'FontWeight', 'bold');
    grid on; axis square; grid minor;
    set(gca, 'FontSize', fs, 'FontWeight', 'bold');

    subplot(1,2,2);
    plot(100*u_spfs, 100*d_bars, 'k-', 'LineWidth' ,lw);
    xlabel('Superficial Gas Velocity (cm/s)', 'FontSize', fs, 'FontWeight', 'bold');
    ylabel('Sauter Mean Diameter (cm)', 'FontSize', fs, 'FontWeight', 'bold');
    grid on; axis square; grid minor;
    set(gca, 'FontSize', fs, 'FontWeight', 'bold');

    %Plot BSD evolution - simple 
    figure();

    subplot(1,2,1);
    plot(dbsf{1}, BSDs_norm{1}, 'b-','LineWidth', lw, 'Color', BSDcmap(1, :)); hold on;
    plot(dbsf{4}, BSDs_norm{4}, 'b-.','LineWidth', lw, 'Color', BSDcmap(4, :))
    plot(dbsf{8}, BSDs_norm{8}, 'b-','LineWidth', lw, 'Color', BSDcmap(8, :))
    plot(dbsf{13}, BSDs_norm{13}, 'b-.','LineWidth', lw, 'Color', BSDcmap(13, :))
    plot(dbsf{18}, BSDs_norm{18}, 'b-','LineWidth', lw, 'Color', BSDcmap(18, :))
    xlabel('Bubble Diameter (m)', 'FontSize', fs);
    ylabel('Bubble Size Distribution (1/m3)', 'FontSize', fs);
    grid on; grid minor; axis square;
    title('Norm. BSD vs. Flow Rate');
    legend('u_{spf} = 0.25 cm/s', 'u_{spf} = 1 cm/s', 'u_{spf} = 5 cm/s', 'u_{spf} = 10 cm/s',...
     'u_{spf} = 15 cm/s', 'FontSize', fs, 'Location', 'northeast');
    set(gca, 'FontSize', fs);

    subplot(1,2,2);
    plot(dbsf{1}, V_PDFs{1}, 'b-','LineWidth', lw, 'Color', BSDcmap(1, :)); hold on;
    plot(dbsf{4}, V_PDFs{4}, 'b-.','LineWidth', lw, 'Color', BSDcmap(4, :))
    plot(dbsf{8}, V_PDFs{8}, 'b-','LineWidth', lw, 'Color', BSDcmap(8, :))
    plot(dbsf{13}, V_PDFs{13}, 'b-.','LineWidth', lw, 'Color', BSDcmap(13, :))
    plot(dbsf{18}, V_PDFs{18}, 'b-','LineWidth', lw, 'Color', BSDcmap(18, :))
    xlabel('Bubble Diameter (m)', 'FontSize', fs);
    ylabel('PDF (mass) (-)', 'FontSize', fs);
    grid on; grid minor; axis square;
    title('Mass PDF vs. Flow Rate');
    legend('u_{spf} = 0.25 cm/s', 'u_{spf} = 1 cm/s', 'u_{spf} = 5 cm/s', 'u_{spf} = 10 cm/s',...
     'u_{spf} = 15 cm/s', 'FontSize', fs, 'Location', 'northeast');
    set(gca, 'FontSize', fs);


%% Test 3 - Full flow rate study - Considering heat transfer and the reaction in liquid tin

%% Plot results from the Full Flow Rate Study

    %Load file
    load('Data/Studies/Flow Rate/flow_rate_sweep_full_2026-04-28-21-45-06.mat');

    %Extract results
    results_all = results_all.results_all;

    %Iterate through cases and analyze
    for i = 1:length(results_all)

        %Pull results and params
        results = results_all{i};
        params = results.params;

        %Extract and plot BSD evolution, gas holdup evolution, final temperature profile, and final conversion profile for each case

    end



end