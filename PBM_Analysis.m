function PBM_Analysis(params)

%% Setup

    %Plot aesthertics
    fs = 18;

    %Load colors
    colors = PlotColors();

%% Comparison of calculated and 

    %Create plot()
    figure(); 
    plot_colors = [colors.trueblue; colors.elecblue;
                   colors.truered; colors.ChiliRed];

    %Setup
    ds = [0.0015, 0.002, 0.003, 0.006];
    iz = 1;
    z = params.cents_y(iz);
    cellinds = ((iz-1).*params.Nms + 1):1:(iz.*params.Nms);
    Ns_cell = ones(size(cellinds));
    uL = 0;

    %Create storage matrix
    BSDs.interp = zeros(length(ds), length(params.fvs_norm_all));
    BSDs.analytic = BSDs.interp;

    %Iterate through size groups
    for id = 1:length(ds)

        %Pull bubble size
        d = ds(id); 

        %Pull epsilon
        turb_eps = 1;
        params.turb.eps(:) = turb_eps;

        %Calculate kolmogorov length scale
        lambda_komogorov = ((params.nus.^3)./turb_eps).^0.25; %m
        lambda_min = 31.4 * lambda_komogorov; %Minimum eddy diameter to break a bubble - integrate up to bubble diameter 

        %Calculate breakage rate and BSD
        [b_eddy, beta, int_ratio] = BreakageEddyAlt(iz, id, d, Ns_cell, lambda_min, params.N_lambdas, params);

        %Calculate using interpolation
        beta_ratio = params.break.beta_ratio{iz}(d, params.turb.eps(iz));
        beta_eddy = params.break.beta{iz}(d, params.turb.eps(iz), params.fvs_norm_all);
        beta_eddy = beta_eddy(:);

        %Store results
        BSDs.analytic(id,:) = beta;
        BSDs.interp(id, :) = beta_eddy;

        %Plot results
        plot(params.fvs_norm_all, BSDs.analytic(id,:), 'LineWidth', 1.5, 'LineStyle', '-', 'Color', plot_colors(id,:)); hold on;
        plot(params.fvs_norm_all, BSDs.interp(id,:), 'LineWidth', 1.5, 'LineStyle', '--', 'Color', plot_colors(id,:));

    end

    %Plot asesthetics
    grid on; grid minor; axis square;
    xlabel('Size (mass or volume) fraction (f_v)'); ylabel('Probability');
    title('Normalized Bubble size distribution (\epsilon = 1 m^2/s^2)');
    ylim([0, max(BSDs.analytic(:)) + 1])
    
    h = zeros(1,6);
    h(1) = plot(NaN, NaN, 'k.', 'MarkerSize', 30, 'Color', plot_colors(1,:));
    h(2) = plot(NaN, NaN, 'k.', 'MarkerSize', 30, 'Color', plot_colors(2,:));
    h(3) = plot(NaN, NaN, 'k.', 'MarkerSize', 30, 'Color', plot_colors(3,:));
    h(4) = plot(NaN, NaN, 'k.', 'MarkerSize', 30, 'Color', plot_colors(4,:));
    h(5) = plot(NaN, NaN, 'k-', 'MarkerSize', 30, 'Color', colors.MITDarkGrey);
    h(6) = plot(NaN, NaN, 'k--', 'MarkerSize', 30, 'Color', colors.MITDarkGrey);

    legend(h, 'd = 1.5 mm', 'd = 2 mm', 'd = 3 mm', 'd = 6 mm', 'Analytic', 'Interp.', 'location', 'north');
    set(gca, 'FontSize', fs);
    %Calculate the same 


    % 
    % %Test cases
    %     eps_ref = 1;
    %     ds_ref = [0.0015, 0.002, 0.003, 0.006];
    %     BSDs = zeros(length(ds_ref), length(params.fvs_norm_all));
    %     for it_debug = 1:length(ds_ref)
    % 
    %         %Calculated interpolated value
    %         BSDs(it_debug,:) = params.betas{iz}(ds_ref(it_debug), eps_ref, params.fvs_norm_all);
    % 
    %         %Renormalize interpolated value
    % 
    % 
    %         %Plot
    %         plot(params.fvs_norm_all, BSDs.analytic(it_debug,:), 'LineWidth', 1.5, 'LineStyle', '-', 'Color', colors(id,:)); hold on;
    %         plot(params.fvs_norm_all, BSDs.interp(id))
    % 
    %     end

%% Single layer convergence example

    %Load colors
    colors = PlotColors();
    fs = 18;

    %Load data files
    data_low = load('Data/Solutions/PBM_output_18-Mar-2026_22-41-27.mat');
    data_high = load('Data/Solutions/PBM_output_18-Mar-2026_22-50-00.mat');

    %Extract final curves
    output_low = data_low.output;
    params_low = output_low.params;
    Y_final_low = output_low.Y(end,:);
    output_high = data_high.output;
    params_high = output_high.params;
    Y_final_high = output_high.Y(end,:);

    %Plot results
    figure();
    subplot(1,2,1);
    plot(params_low.dms, Y_final_low, 'LineWidth', 1.5, 'LineStyle', '-', 'Color', colors.trueblue); hold on;
    plot(params_high.dms, Y_final_high, 'ko', 'LineWidth', 1.5);
    grid on; grid minor; axis square;
    %xlabel('Bubble diameter (m)'); ylabel('Numeric Density (1/m^3)');
    title('Single Layer Convergence Example');
    ylim([0, 18E4]);
    xlabel('Diameter (m)'); ylabel('Numeric Density')
    legend('Di = 6 mm', 'Di = 20 mm');
    set(gca, 'FontSize', fs, 'FontWeight', 'bold');


%% Identify the best empirical coefficient for 

    %Define cases
    NKs = 15;
    Kws = logspace(log10(0.01), log10(20),NKs);

    %Create storage vector for outstructs
    outstructs = cell(1,length(Kws));

    %Specify output folder
        output_folder = 'Data/Studies/Coeff_Analysis/';


    parfor ik = 1:NKs
         
        %Define water structure
        water = struct();
        water.dyn_visc = @(T) 0.0010016 .* ones(size(T)); %Pa*s - at 20 C
        water.mol_mass = 0.0180153; %kg/mol
        water.surf_tension =  @(T) 0.0728 .*ones(size(T)); 
        water.density = @(T) 998.21 .* ones(size(T)); 
        water.Cp = @(T) 4157 .* ones(size(T)); 
        water.therm_cond = @(T) 0.598 .* ones(size(T)); 
        water.name = 'water';

        %Test conditions
        u_spfs = [0.01, 0.03, 0.05, 0.12, 0.16];
        alpha_gs = [0.038 0.109, 0.154, 0.198, 0.225];
        epss = 10 .* u_spfs;
    
        %Create storage array
        results = cell(1, length(u_spfs));
        ys_out = cell(1, length(u_spfs));
    
        ys_final = cell(1,length(u_spfs));

        outstruct = struct();
        
    
        %Iterate through cases
        for ic = 1:length(u_spfs)
    
            %Initialize inputs and specify coalescence/breakage kernel
            inputs = PBM_inputs();
            inputs.src.coalesce_model = 'Wang_2005';
            inputs.src.breakage_model = 'Wang_2005';
    
            %Modify inputs
            
            inputs.disc.Nms = 40;
            inputs.disc.mesh_hybrid_cells = round(inputs.disc.Nms .* 0.5);  
            inputs.reactor.liquid = water;
            inputs.reactor.u_spf_orifice = u_spfs(ic); %m/s
            inputs.reactor.T_gas_i_mu = 20 + 273.15;
            inputs.reactor.T = 20 + 273.15;
            inputs.reactor.T_orifice = inputs.reactor.T;
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
            inputs.src.c_wake = Kws(ik);
            inputs.sim.t_end = 1;
            inputs.sol.src_delay = 0;
            inputs.sol.break_file = 'break_water_Nd-20_Nz-1_Ne-20_TL-293.mat';
            inputs.sol.postprocess = false; %Disabled to prevent figures 
    
            %Run PBM
            results{ic} = PBM_v3(inputs);
            ys_out{ic} = results{ic}.y;
            ys_final{ic} = results{ic}.y(end,:);
            
        end
    
        %Log results
        outstruct.results = results; outstruct.ys_out = ys_out;
        outstruct.ys_final = ys_final;
        outstruct.Kws = Kws;
        outstructs{ik} = outstruct;


    end

    %Log other values
    


    %Save data
    
    timestamp = datestr(now, 'yyyy-mm-dd-HH-MM-SS');
    outfilename = sprintf('flowrate_coeffanalysis_%s.mat', timestamp);
    outfilepath = fullfile(output_folder, outfilename);
    save(outfilepath, 'outstructs');

%% Plot the results from the coeff study 

    %Load data
    data = load('Data/Studies/Coeff_Analysis/flowrate_coeffanalysis_2026-04-24-18-35-30.mat');
    outstructs = data.outstructs;
    Kws = outstructs{1}.Kws;

    %
    refcolors = PlotColors();

    %Specify plot colors
    colors = {[0.85, 0.3250, 0.0980];
              [0.929, 0.694, 0.125];
             [0.494, 0.184, 0.556];
              [0.466, 0.674, 0.188];
               [0.301, 0.745, 0.933]};

    %Create overarching figure for comparing results across cases
    kfig = figure();

    %Iterate through cases
    for ic = 1:length(outstructs)

        %Pull results
        outstruct = outstructs{ic};
        results = outstruct.results;

        %Create figure for pltos
        cfig = figure();

        %Iterate through results
        for ir = 1:length(results)

            %Pull result
            result = results{ir};
            

            %Isolate relevant values
            params = result.params;
            y_final = result.y(end,2:end); 

            %Calculate equivalent diameters
            p = params.p_surf;
            T_bar = params.T_liq;
            ns_gas = params.nms;
            Vs = (ns_gas * params.R * T_bar)/p;
            ds = (6.*Vs/pi).^(1/3);
            ds = ds(2:end);
            mms = params.mms(2:end);


            %Calculate volume PDF 
            Vbs = Vs(2:end) .* y_final;
            V_total = sum(Vbs);
            V_PDF = Vbs./V_total;
            V_PDF_sum = sum(V_PDF);
            V_PDF_trapz = trapz(1000*ds, V_PDF);
            V_PDF = V_PDF/V_PDF_trapz;

            %Plot all profiles for this case
            figure(cfig)
            subplot(2,1,1);
            plot(ds, y_final, 'LineWidth', 2, 'Color', colors{ir}); hold on;

            subplot(2,1,2);
            plot(ds, V_PDF, 'LineWidth', 2, 'Color', colors{ir}); hold on;

            %Plot for one case
            if ir == 3

                %Change figures
                figure(kfig);

                if ic == 1
                    y_final_init = y_final;
                    V_PDF_init = V_PDF;
                elseif ic == length(results)
                    plotcolor = refcolors.hydrogen;
                    lw = 3;
                    subplot(1,2,1);
                    plot(ds, y_final_init, 'LineWidth', lw, 'Color', plotcolor); hold on;

                    subplot(1,2,2);
                    plot(ds, V_PDF_init, 'LineWidth', lw, 'Color', plotcolor); hold on;  

                end

                %Change figures
                figure(kfig);

                %Select color
                if ic == 1
                    plotcolor = refcolors.hydrogen;
                    lw = 3;
                elseif ic == length(outstructs)
                    plotcolor = refcolors.ChiliRed;
                    lw = 3;
                else
                    plotcolor = refcolors.MITDarkGrey;
                    lw = 1;
                end

                %Plot
                subplot(1,2,1);
                plot(ds, y_final, 'LineWidth', lw, 'Color', plotcolor); hold on;

                subplot(1,2,2);
                plot(ds, V_PDF, 'LineWidth', lw, 'Color', plotcolor); hold on;  

                
            end

        end

        %Label the plots
        figure(cfig);
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



    %Plot aesthetics for
    figure(kfig);
    subplot(1,2,1);
    xlabel('Bubble Diameter (m)'); ylabel('Numeric Density (1/m3)');
    grid on; grid minor; axis square; 
    set(gca, 'FontSize', 14, 'FontWeight', 'bold');
    title('u = 0.05 m/s', 'FontSize', 16, 'FontWeight', 'bold');

    %Create legend
    h = zeros(1, 3);
    h(1) = plot(NaN, NaN, 'k-', 'LineWidth', 2, 'Color', refcolors.hydrogen);
    h(2) = plot(NaN, NaN, 'k-', 'LineWidth', 2, 'Color', refcolors.MITDarkGrey);
    h(3) = plot(NaN, NaN, 'k-', 'LineWidth', 2, 'Color', refcolors.ChiliRed);
    legend(h, sprintf('K = %0.4f', Kws(1)), 'Other K', sprintf('K = %0.4f', Kws(end)), 'location', 'northeast');

    figure(kfig);
    subplot(1,2,2);
    xlabel('Bubble Diameter (m)'); ylabel('PDF (volume) of d_b (1/m)');
    grid on; grid minor; axis square; 
    set(gca, 'FontSize', 14, 'FontWeight', 'bold');
    title('u = 0.05 m/s', 'FontSize', 16, 'FontWeight', 'bold');

    %Create legend
    h = zeros(1, 3);
    h(1) = plot(NaN, NaN, 'k-', 'LineWidth', 2, 'Color', refcolors.hydrogen);
    h(2) = plot(NaN, NaN, 'k-', 'LineWidth', 2, 'Color', refcolors.MITDarkGrey);
    h(3) = plot(NaN, NaN, 'k-', 'LineWidth', 2, 'Color', refcolors.ChiliRed);
    legend(h, sprintf('K = %0.4f', Kws(1)), 'Other K', sprintf('K = %0.4f', Kws(end)), 'location', 'northeast');


end