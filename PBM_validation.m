function [] = PBM_validation()

%% Setup


    %Load material properties
    addpath('Material Properties/');
    load tin.mat
    load hydrogen.mat
    load carbon.mat
    load methane.mat
    load argon.mat
    load water.mat


%% Besagni et al. (2019) cases - run

    %Specify test cases
    ugs = [0.0037, 0.0111, 0.0188]; %[0.0037, 0.0074, 0.0111, 0.0144, 0.0188]; %Superficial gas velocities (m/s)

    %Initialize inputs
    inputs = PBM_inputs();
    inputs.sim.t_end = 10;
    inputs.inlet.d.mu_i = 0.002;
    inputs.inlet.d.std_i = 0.0005;
    inputs.disc.Nms = 60;
    inputs.disc.r_max = 0.01; %m - 
    inputs.disc.r_min = 0.0005;
    inputs.sol.heat.active = false;
    inputs.sol.single_layer = true;
    inputs.reactor.liquid = water;
    inputs.reactor.T = 298.15; %1200 + 273.15;
    inputs.reactor.T_orifice = inputs.reactor.T;
    inputs.reactor.T_gas_i_mu = 293.15;
    

    %Allocate results structure
    results = struct([]);

    %Run cases 
    for ig = 1:length(ugs)

        %Update inputs
        inputs.reactor.u_spf_orifice = ugs(ig);
        results(ig).u_spf = ugs(ig);

        % %Run with Wang models
        inputs.src.coalesce_model = 'Wang_2005';
        inputs.src.breakage_model = 'Wang_2005';
        results(ig).wang_wang = PBM_v3(inputs);

        % %Run with Wang coalesce, luo breakage
        inputs.src.coalesce_model = 'Wang_2005';
        inputs.src.breakage_model = 'Luo_Svendson_1996';
        results(ig).wang_luo = PBM_v3(inputs);
        
        % %Run with Prince coalesce, Luo breakage
        inputs.src.coalesce_model = 'PrinceBlanch_1990';
        inputs.src.breakage_model = 'Luo_Svendson_1996';
        results(ig).prince_luo = PBM_v3(inputs);
         
        %Run with Prince coalesce, Wang breakage
        inputs.src.coalesce_model = 'PrinceBlanch_1990';
        inputs.src.breakage_model = 'Wang_2005';
        results(ig).prince_wang = PBM_v3(inputs);


        % %Store results 
        % if ig == 1
        %     results_0037.wang_2005 = results;
        % elseif ig == 2
        %     results_0074 = results;
        % elseif ig == 3
        %     results_0111 = results;
        % end

        
        %Compare to data
        x = 1

    end

    %Define output 
    %results.results_0037 = results_0037;
    %results.results_0074 = results_0074;
    %results.results_0111 = results_0111;

    %Define output filename and path, using a current timestamp to avoid overwriting previous results
    timestamp = datestr(now, 'yyyy-mm-dd_HH-MM');
    outname = sprintf('Besagni_2019_sim_results_Nm-%d_%s.mat',inputs.disc.Nms, timestamp);
    %outname = [ timestamp, '.mat'];
    outpath = ['Validation/Besagni_2019/', outname];

    %Save outpuit
    save(outpath, 'results');


%% Besagni et al. (2019) cases - plot 

    %Define plot aesthetics
    lw = 2; %Line width
    fs = 18; %Font size
    ms = 18;

    %Load results
    if ~exist('outname')
        results = load('Validation/Besagni_2019/Besagni_2019_sim_results_Nm-60_v4.mat');
        results = results.results;
    end

    %Iterate through u_spfs and plot results for each case
    N_uspfs = length(results);
    for iu = 1:N_uspfs

        %Extract gas holdup
        u_spf = results(iu).u_spf;

        %Define 
        if u_spf == 0.0037
            uspftxt = '0037';
        elseif u_spf == 0.0074
            uspftxt = '0074';
        elseif u_spf == 0.0111
            uspftxt = '0111';
        elseif u_spf == 0.0149
            uspftxt = '0144';
        elseif u_spf == 0.0188            
            uspftxt = '0188';
        end

        %Load validation data for this case
        filename = sprintf('Validation/Besagni_2019/Besagni_2019_ug_%s.csv', uspftxt);
        valdata = ExtractExpVals(readtable(filename));

        %Create figure for this case
        figname = sprintf('ug = %.4f m/s', u_spf);
        figure('Name', figname, 'units', 'normalized', 'outerposition', [0 0 10/16 1]);

        %Plot Wang wang case
        subplot(2,2,1);
        wang_norm = results(iu).wang_wang.y(end,:)./trapz(1000*results(iu).wang_wang.dbsout, results(iu).wang_wang.y(end,:));   
        plot(results(iu).wang_wang.dbsout.*1000, wang_norm, 'b-', 'LineWidth', lw); hold on;
        plot(valdata.xs, valdata.ys, 'r.', 'MarkerSize', ms);
        xlabel('Diameter (mm)'); ylabel('Number Density (1/m^3)'); title(sprintf('Wang coalescence, Wang breakage'));
        grid on; grid minor;    axis square;
        legend('PBM', 'Experiments');

        %Plot Wang luo case
        subplot(2,2,2);
        wang_luo_norm = results(iu).wang_luo.y(end,:)./trapz(1000*results(iu).wang_luo.dbsout, results(iu).wang_luo.y(end,:));   
        plot(results(iu).wang_luo.dbsout.*1000, wang_luo_norm, 'b-', 'LineWidth', lw); hold on;
        plot(valdata.xs, valdata.ys, 'r.', 'MarkerSize', ms);
        xlabel('Diameter (mm)'); ylabel('Number Density (1/m^3)'); title(sprintf('Wang coalescence, Luo breakage'));
        grid on; grid minor;    axis square;
        legend('PBM', 'Experiments');   

        %Plot Prince luo case
        subplot(2,2,3);
        prince_luo_norm = results(iu).prince_luo.y(end,:)./trapz(1000*results(iu).prince_luo.dbsout, results(iu).prince_luo.y(end,:));   
        plot(results(iu).prince_luo.dbsout.*1000, prince_luo_norm, 'b-', 'LineWidth', lw); hold on;
        plot(valdata.xs, valdata.ys, 'r.', 'MarkerSize', ms);
        xlabel('Diameter (mm)'); ylabel('Number Density (1/m^3)'); title(sprintf('Prince coalescence, Luo breakage'));
        grid on; grid minor;    axis square;
        legend('PBM', 'Experiments');       

        %Plot Prince wang case
        subplot(2,2,4);
        prince_wang_norm = results(iu).prince_wang.y(end,:)./trapz(1000*results(iu).prince_wang.dbsout, results(iu).prince_wang.y(end,:));   
        plot(results(iu).prince_wang.dbsout.*1000, prince_wang_norm, 'b-', 'LineWidth', lw); hold on;
        plot(valdata.xs, valdata.ys, 'r.', 'MarkerSize', ms);
        xlabel('Diameter (mm)'); ylabel('Number Density (1/m^3)'); title(sprintf('Prince coalescence, Wang breakage'));
        grid on; grid minor;    axis square;
        legend('PBM', 'Experiments');

    end






end


function results = ExtractExpVals(results_in)

    %Pull values
    xs_in =results_in(:,1); xs_out = table2array(xs_in);
    ys_in = results_in(:,2); ys_out = table2array(ys_in);

    % %Identify points where slope changes
    % dxs = diff(xs_in);
    % dys = diff(ys_in);
    % ms = dys./dxs; %Calculate slope between each point


    %Define output
    results.xs = xs_out;
    results.ys = ys_out;
    

end