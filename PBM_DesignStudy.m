% This function studys the design space for the bubble column to 
% identify a narrower range for later optimization

function PBM_DesignStudy()

%% Setup

    %Settings
    debug = true;
    saveplots = true;

    %Plot settings
    lw = 2; %Line width
    ms = 8; %Marker size
    fs = 14; %Font size

    %Material Properties
    addpath('Material Properties/');
    load tin.mat
    load hydrogen.mat
    load carbon.mat
    load methane.mat
    load argon.mat

    %Update values
    tin.dyn_visc = @(T) 0.00031 .* exp(6.171./(0.008314 .* T));
    tin.name = "tin";



%% Gas Flow Rate and Temperature Study - Analysis 

    %Specify output folder
    output_folder = 'Data/Studies/Flow Rate/';

    %Test conditions
    Ts = 1100:100:1300;
    u_spfs = [0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1];
   



    %Create storage array
    results = cell(length(Ts), length(u_spfs));
    ys_out = cell(length(Ts), length(u_spfs));
    ys_final = cell(length(Ts), length(u_spfs));
    %Ts_in
    

    

    %Iterate through cases
    parfor ic = 1:length(u_spfs)
        for it = 1:3
    
            %Initialize inputs and specify coalescence/breakage kernel
            inputs = PBM_inputs();
    
            %Modify inputs - simulation
            inputs.sim.t_end = 3;
    
            %Modify inpouts - source terms
            inputs.src.coalesce_model = 'Wang_2005';
            inputs.src.breakage_model = 'Wang_2005';
            inputs.src.c_wake = 1.2;
            inputs.sol.src_delay = 0;
            inputs.src.breakage_active = true;
            inputs.src.coalesce_active = true;
    
            %Modify
    
            %Modify inputs - discretization and solver
            inputs.disc.Nms = 50;
            inputs.disc.Nz = 8;
            inputs.disc.mesh_hybrid_cells = round(inputs.disc.Nms .* 0.5);  
            inputs.sol.heat.active = true;
            inputs.sol.react.active = true;
            inputs.sol.solve_ub = true;
            inputs.sol.single_layer = false;
            inputs.sol.postprocess = false; %Disabled to prevent figures 
            
    
            %Modify inputs - geometry and physical conditions
            inputs.reactor.H = 0.2;
            inputs.reactor.liquid = tin;   
            inputs.reactor.u_spf_orifice = u_spfs(ic); %m/s
            inputs.reactor.T_gas_i_mu = 1073.15;
            inputs.reactor.T = Ts(it) + 273.15;
            inputs.reactor.T_orifice = inputs.reactor.T;
            inputs.reactor.T_min_i = 700 + 273.15;
            inputs.inlet.d.mu_i = 0.01;
            inputs.inlet.d.std_i = 0.001;
            
    
            %Run PBM
            output = PBM_v3(inputs);

            %Store results
            results{it,ic} = output;
            ys_out{it,ic} = output.y;
            ys_final{it,ic} = output.y(end,:);
            
        end

    end

    %Save data
    outstruct.results = results; outstruct.ys_out = ys_out;
    outstruct.ys_final = ys_final;
    timestamp = datestr(now, 'yyyy-mm-dd-HH-MM-SS');
    outfilename = sprintf('design_%s.mat', timestamp);
    outfilepath = fullfile(output_folder, outfilename);
    save(outfilepath, 'outstruct');


%% Gas Flow Rate and Temperature Study - Plots

    %Load file
    results_all = load('Data/Studies/Flow Rate/design_2026-05-13-08-06-43.mat');

    %Extract data
    results = results_all.outstruct.results;
    ys_out = results_all.outstruct.ys_out;
    ys_final = results_all.outstruct.ys_final;

    %Plot settings 
    lw = 2;
    fs = 16;
    plotcolors = PlotColors();

    %Pull Case values - update in future 
    Ts = 1100:100:1300;
    u_spfs = [0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.125, 0.15];

    %Create color map
    BSDcmap = nebula(length(u_spfs));

    %Allocate storage vectors
    alphags = zeros(length(Ts), length(u_spfs));
    Tbars = zeros(length(Ts), length(u_spfs));
    Xbars = zeros(length(Ts), length(u_spfs));
    dbars = zeros(length(Ts), length(u_spfs));
    BSDs = cell(length(Ts), length(u_spfs));
    BSDs_norm = cell(length(Ts), length(u_spfs));
    V_PDFs = cell(length(Ts), length(u_spfs));
    dbsf = cell(length(Ts), length(u_spfs));

    %Iterate through results and pull values
    for it = 1:length(Ts)
        for iu = 1:length(u_spfs)

            %Pull data
            output = results{it,iu};
            y_out = ys_out{it,iu};
            y_final = ys_final{it,iu};
            params = output.params;

            %Pull values
            alphags(it,iu) = output.alphag_avg;
            Tbars(it,iu) = output.Tg_mix(end);
            Xbars(it,iu) = output.Xg_mix(end);
            dbars(it,iu) = output.dsauter(end);
            dbsf{it, iu} = output.dbsout(end,:);
            BSDs{it,iu} = output.mdists(end,:);
            BSDs_norm{it,iu} = BSDs{it,iu}./trapz(dbsf{it,iu}, BSDs{it,iu});

            %Calculate mass/volume probability distribution
            Vbs = params.Vms(1:end) .* BSDs{it,iu};
            V_total = sum(Vbs);
            V_PDF = Vbs./V_total;
            V_PDF_sum = sum(V_PDF);
            V_PDF_trapz = trapz(1000*dbsf{it,iu}, V_PDF);
            V_PDF = V_PDF/V_PDF_trapz;
            V_PDFs{it,iu} = V_PDF;

        end

        %Plot BSDs and PDFs
        figure();

        subplot(1,2,1);
        plot(dbsf{it,2}, BSDs_norm{it,2}, 'b-','LineWidth', lw, 'Color', BSDcmap(1, :)); hold on;
        plot(dbsf{it,6}, BSDs_norm{it,6}, 'b-.','LineWidth', lw, 'Color', BSDcmap(4, :))
        plot(dbsf{it,8}, BSDs_norm{it,8}, 'b-','LineWidth', lw, 'Color', BSDcmap(7, :))
        plot(dbsf{it,10}, BSDs_norm{it,10}, 'b-.','LineWidth', lw, 'Color', BSDcmap(10, :))
        xlabel('Bubble Diameter (m)', 'FontSize', fs);
        ylabel('Bubble Size Distribution (1/m3)', 'FontSize', fs);
        grid on; grid minor; axis square;
        title('Norm. BSD vs. Flow Rate');
        legend('u_{spf} = 1 cm/s', 'u_{spf} = 5 cm/s', 'u_{spf} = 10 cm/s',...
         'u_{spf} = 15 cm/s', 'FontSize', fs, 'Location', 'northeast');
    
        subplot(1,2,2);
        plot(dbsf{it,2}, V_PDFs{it,2}, 'b-','LineWidth', lw, 'Color', BSDcmap(1, :)); hold on;
        plot(dbsf{it,6}, V_PDFs{it,6}, 'b-.','LineWidth', lw, 'Color', BSDcmap(4, :))
        plot(dbsf{it,8}, V_PDFs{it,8}, 'b-','LineWidth', lw, 'Color', BSDcmap(7, :))
        plot(dbsf{it,10}, V_PDFs{it,10}, 'b-.','LineWidth', lw, 'Color', BSDcmap(10, :))

        xlabel('Bubble Diameter (m)', 'FontSize', fs);
        ylabel('PDF (mass) (-)', 'FontSize', fs);
        grid on; grid minor; axis square;
        title('Mass PDF vs. Flow Rate');
        legend('u_{spf} = 1 cm/s', 'u_{spf} = 5 cm/s', 'u_{spf} = 10 cm/s',...
         'u_{spf} = 15 cm/s', 'FontSize', fs, 'Location', 'northeast');


    end

    %Specify color progression

    %Plot temperature and conversion versus superficial gas velocity
    figure();

    subplot(2,2,1);


    yyaxis left;
    plot(100.*u_spfs, Tbars(1,:)-273.15, 'b--','LineWidth', lw, 'Color', plotcolors.ChiliRed); hold on;
    ylabel('Gas Exit Temp. (\circ C)', 'FontSize', fs);
    set(gca,'YColor', plotcolors.ChiliRed);

    yyaxis right;
    plot(100.*u_spfs, 100.*Xbars(1,:), 'b--','LineWidth', lw, 'Color', plotcolors.hydrogen); hold on;
    ylabel('Conversion (%)', 'FontSize', fs);
    set(gca,'YColor', plotcolors.hydrogen);
    title('T_{liquid} = 1100 \circ C Case')
    xlabel('Gas Velocity (cm/s)', 'FontSize', fs);
    grid on; grid minor; 
    xlim([0,10]);

    subplot(2,2,2);
    yyaxis left;
    plot(100.*u_spfs, Tbars(2,:)-273.15, 'b-.','LineWidth', lw, 'Color', plotcolors.ChiliRed); hold on;
    ylabel('Gas Exit Temp. (\circ C)', 'FontSize', fs);
    set(gca,'YColor', plotcolors.ChiliRed);
    
    yyaxis right;
    plot(100.*u_spfs, 100.*Xbars(2,:), 'b-.','LineWidth', lw, 'Color', plotcolors.hydrogen); hold on;
    ylabel('Conversion (%)', 'FontSize', fs);
    set(gca,'YColor', plotcolors.hydrogen);
    title('T_{liquid} = 1200 \circ C Case')
    xlabel('Gas Velocity (cm/s)', 'FontSize', fs);
    grid on; grid minor; 
    xlim([0,10]);

    subplot(2,2,[3]);
    yyaxis left;
    plot(100.*u_spfs, Tbars(3,:)-273.15, 'b-','LineWidth', lw, 'Color', plotcolors.ChiliRed); hold on;
    ylabel('Gas Exit Temp. (\circ C)', 'FontSize', fs);
    set(gca,'YColor', plotcolors.ChiliRed);

    yyaxis right;

    plot(100.*u_spfs, 100.*Xbars(3,:), 'b-','LineWidth', lw, 'Color', plotcolors.hydrogen);
    ylabel('Conversion (%)', 'FontSize', fs);
    set(gca,'YColor', plotcolors.hydrogen);
    title('T_{liquid} = 1300 \circ C Case')
    xlabel('Gas Velocity (cm/s)', 'FontSize', fs);
    grid on; grid minor; 
    xlim([0,10]);

    subplot(2,2,4);
    plot(100.*u_spfs, 100.*Xbars(1,:), 'b--','LineWidth', lw, 'Color', plotcolors.hydrogen); hold on;
    plot(100.*u_spfs, 100.*Xbars(2,:), 'b-.','LineWidth', lw, 'Color', plotcolors.hydrogen); hold on;
    plot(100.*u_spfs, 100.*Xbars(3,:), 'b-','LineWidth', lw, 'Color', plotcolors.hydrogen);
    ylabel('Conversion (%)', 'FontSize', fs);
    title('Conversion (%) - All Cases')
    xlabel('Gas Velocity (cm/s)', 'FontSize', fs);
    grid on; grid minor; 
    ylim([0, 105]);
    xlim([0,10]);

    %Plot gas holdup and bubble diameter versus superficial gas velocity
    figure();

    subplot(1,2,1);
    plot(100.*u_spfs, 100.*alphags(1,:), 'k-','LineWidth', lw, 'Color', plotcolors.hydrogen); hold on;
    plot(100.*u_spfs, 100.*alphags(2,:), 'k-','LineWidth', lw, 'Color', plotcolors.MITDarkGrey); hold on;
    plot(100.*u_spfs, 100.*alphags(3,:), 'k-','LineWidth', lw, 'Color', plotcolors.ChiliRed); hold on;
    grid on; grid minor; axis square;
    xlabel('Gas Velocity (cm/s)', 'FontSize', fs);
    title('Gas Holdup vs. Superficial Velocity');
    ylabel('Gas Holdup (%)', 'FontSize', fs);
    xlim([0,10]);
    legend('T_{liquid} = 1100 \circC', 'T_{liquid} = 1200 \circC', 'T_{liquid} = 1300 \circC', 'FontSize', fs, 'Location', 'northeast');

    subplot(1,2,2);
    plot(100.*u_spfs, 1000.*dbars(1,:), 'b-','LineWidth', lw, 'Color', plotcolors.hydrogen); hold on;
    plot(100.*u_spfs, 1000.*dbars(2,:), 'b-','LineWidth', lw, 'Color', plotcolors.MITDarkGrey);
    plot(100.*u_spfs, 1000.*dbars(3,:), 'b-','LineWidth', lw, 'Color', plotcolors.ChiliRed);
    grid on; grid minor; axis square;
    xlabel('Gas Velocity (cm/s)', 'FontSize', fs);
    ylabel('Sauter Mean Diameter (mm)', 'FontSize', fs);  
    title('Mean Diameter vs. Superficial Velocity');
    legend('T_{liquid} = 1100 \circC', 'T_{liquid} = 1200 \circC', 'T_{liquid} = 1300 \circC', 'FontSize', fs, 'Location', 'northeast');
    xlim([0,10]);

end