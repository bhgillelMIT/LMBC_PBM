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

    water.name = 'water';


%% Besagni et al. (2019) cases - run

    %Specify test cases
    ugs = [0.0037, 0.0111, 0.0188]; %[0.0037, 0.0074, 0.0111, 0.0144, 0.0188]; %Superficial gas velocities (m/s)
    alpha_gs = [0.039, 0.085, 0.129];
    epss = [0.092, 0.246, 0.421];

    %Specify other static inputs
    Nms = 60;

    % %Initialize inputs
    % inputs = PBM_inputs();
    % inputs.sim.t_end = 10;
    % inputs.inlet.d.mu_i = 0.002;
    % inputs.inlet.d.std_i = 0.0005;
    % inputs.disc.Nms = 60;
    % inputs.disc.r_max = 0.01; %m - 
    % inputs.disc.r_min = 0.0005;
    % inputs.sol.heat.active = false;
    % inputs.sol.single_layer = true;
    % inputs.reactor.liquid = water;
    % inputs.reactor.T = 298.15; %1200 + 273.15;
    % inputs.reactor.T_orifice = inputs.reactor.T;
    % inputs.reactor.T_gas_i_mu = 293.15;
    

    %Allocate results structure
    results = struct([]);

    %Run cases 
    parfor ig = 1:length(ugs)

        %Initialize inputs
        inputs = PBM_inputs();
        inputs.sim.t_end = 10;
        inputs.inlet.d.mu_i = 0.002;
        inputs.inlet.d.std_i = 0.0005;
        inputs.disc.Nms = Nms;
        inputs.disc.r_max = 0.01; %m - 
        inputs.disc.r_min = 0.00001;
        inputs.sol.heat.active = false;
        inputs.sol.react.active = false;
        inputs.sol.single_layer = true;
        inputs.reactor.liquid = water;
        inputs.reactor.T = 298.15; %1200 + 273.15;
        inputs.reactor.T_orifice = inputs.reactor.T;
        inputs.reactor.T_gas_i_mu = 293.15;

        if true
            inputs.src.eps_manual = epss(ig);
            inputs.src.solve_eps = false;
            inputs.src.alphag_manual = alpha_gs(ig);
            inputs.src.solve_alphag = false;
        end

        %Update inputs for case
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
    outname = sprintf('Besagni_2019_sim_results_Nm-%d_%s.mat',Nms, timestamp);
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
        results = load('Validation/Besagni_2019/Besagni_2019_sim_results_Nm-60_v5.mat');
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
        xlabel('Diameter (mm)'); ylabel('PDF Number (-)'); title(sprintf('Wang coalescence, Wang breakage'));
        grid on; grid minor;    axis square;
        legend('PBM', 'Experiments');

        %Plot Wang luo case
        subplot(2,2,2);
        wang_luo_norm = results(iu).wang_luo.y(end,:)./trapz(1000*results(iu).wang_luo.dbsout, results(iu).wang_luo.y(end,:));   
        wang_luo_norm(1) = 0;
        plot(results(iu).wang_luo.dbsout.*1000, wang_luo_norm, 'b-', 'LineWidth', lw); hold on;
        plot(valdata.xs, valdata.ys, 'r.', 'MarkerSize', ms);
        xlabel('Diameter (mm)'); ylabel('PDF Number (-)'); title(sprintf('Wang coalescence, Luo breakage'));
        grid on; grid minor;    axis square;
        legend('PBM', 'Experiments');   

        %Plot Prince luo case
        subplot(2,2,3);
        prince_luo_norm = results(iu).prince_luo.y(end,:)./trapz(1000*results(iu).prince_luo.dbsout, results(iu).prince_luo.y(end,:));
        prince_luo_norm(1) = 0;
        plot(results(iu).prince_luo.dbsout.*1000, prince_luo_norm, 'b-', 'LineWidth', lw); hold on;
        plot(valdata.xs, valdata.ys, 'r.', 'MarkerSize', ms);
        xlabel('Diameter (mm)'); ylabel('PDF Number (-)'); title(sprintf('Prince coalescence, Luo breakage'));
        grid on; grid minor;    axis square;
        legend('PBM', 'Experiments');       

        %Plot Prince wang case
        subplot(2,2,4);
        prince_wang_norm = results(iu).prince_wang.y(end,:)./trapz(1000*results(iu).prince_wang.dbsout, results(iu).prince_wang.y(end,:));   
        plot(results(iu).prince_wang.dbsout.*1000, prince_wang_norm, 'b-', 'LineWidth', lw); hold on;
        plot(valdata.xs, valdata.ys, 'r.', 'MarkerSize', ms);
        xlabel('Diameter (mm)'); ylabel('PDF Number (-)'); title(sprintf('Prince coalescence, Wang breakage'));
        grid on; grid minor;    axis square;
        legend('PBM', 'Experiments');

    end


%% Wang et al (2005) plots

    %Define plot aesthetics
    lw = 1; %Line width
    fs = 18; %Font size
    ms = 8;

    %Test cases 
    ugs = [0.016, 0.045, 0.077]; %[0.0037, 0.0074, 0.0111, 0.0144, 0.0188]; %Superficial gas velocities (m/s)
    alpha_gs = [0.039, 0.085, 0.129];
    epss = [0.092, 0.246, 0.421];

    %Load data 
    pbmdata = load('Validation/Wang_2005/Wang_2005_sim_results_Nm-60_v1.mat');
    data1 = readtable('Validation/Wang_2005/Wang_2005_ug0016_ul017.csv');
    data2 = readtable('Validation/Wang_2005/Wang_2005_ug0045_ul023.csv');
    data3 = readtable('Validation/Wang_2005/Wang_2005_ug0077_ul025.csv');

    results = pbmdata.results;

    %Create figure
    figure('Name', 'Wang Comparison', 'units', 'normalized', 'outerposition', [0 0 10/16 1]);

    %Iterate through cases 
    N_uspfs = length(results);
    for iu = 1:N_uspfs

        %Extract gas holdup
        u_spf = ugs(iu);

        %Define 
        if iu == 1
            refdata = data1;
            mktype = "square";
            lntype = '-';
            clr = 'k';
        elseif iu == 2
            refdata = data2;
            lntype = ':';
            mktype = 'o';
            clr = 'k';
        elseif iu == 3
            refdata = data3;
            lntype = '-.';
            mktype = '^';
            clr = 'k';
        end


        lnspec = [clr, lntype];

        %Create figure for this case
      %  figname = sprintf('ug = %.4f m/s', u_spf);
       % figure('Name', figname, 'units', 'normalized', 'outerposition', [0 0 10/16 1]);


        valdata = ExtractExpVals(refdata);

        %Plot Wang wang case
        subplot(2,2,1);
        wang_norm = results(iu).wang_wang.y(end,:)./trapz(1000*results(iu).wang_wang.dbsout, results(iu).wang_wang.y(end,:));   
        plot(results(iu).wang_wang.dbsout.*1000, wang_norm, lnspec, 'LineWidth', lw); hold on;
        plot(valdata.xs, valdata.ys, 'k.', 'MarkerSize', ms, 'Marker', mktype);
        xlabel('Diameter (mm)'); ylabel('PDF Number (-)'); title(sprintf('Wang coalescence, Wang breakage'));
        grid on; grid minor;    axis square;
        legend('PBM', 'Experiments');
        xlim([0,10])

        %Plot Wang luo case
        subplot(2,2,2);
        wang_luo_norm = results(iu).wang_luo.y(end,:)./trapz(1000*results(iu).wang_luo.dbsout, results(iu).wang_luo.y(end,:));   
        wang_luo_norm(1) = 0;
        plot(results(iu).wang_luo.dbsout.*1000, wang_luo_norm, lnspec, 'LineWidth', lw); hold on;
        plot(valdata.xs, valdata.ys, 'k.', 'MarkerSize', ms, 'Marker', mktype);
        xlabel('Diameter (mm)'); ylabel('PDF Number (-)'); title(sprintf('Wang coalescence, Luo breakage'));
        grid on; grid minor;    axis square;
        legend('PBM', 'Experiments');  
        xlim([0,10])

        %Plot Prince luo case
        subplot(2,2,3);
        prince_luo_norm = results(iu).prince_luo.y(end,:)./trapz(1000*results(iu).prince_luo.dbsout, results(iu).prince_luo.y(end,:));
        prince_luo_norm(1) = 0;
        plot(results(iu).prince_luo.dbsout.*1000, prince_luo_norm, lnspec, 'LineWidth', lw); hold on;
        plot(valdata.xs, valdata.ys, 'k.', 'MarkerSize', ms, 'Marker', mktype);
        xlabel('Diameter (mm)'); ylabel('PDF Number (-)'); title(sprintf('Prince coalescence, Luo breakage'));
        grid on; grid minor;    axis square;
        legend('PBM', 'Experiments');   
        xlim([0,10])

        %Plot Prince wang case
        subplot(2,2,4);
        prince_wang_norm = results(iu).prince_wang.y(end,:)./trapz(1000*results(iu).prince_wang.dbsout, results(iu).prince_wang.y(end,:));   
        plot(results(iu).prince_wang.dbsout.*1000, prince_wang_norm, lnspec, 'LineWidth', lw); hold on;
        plot(valdata.xs, valdata.ys, 'k.', 'MarkerSize', ms, 'Marker', mktype);
        xlabel('Diameter (mm)'); ylabel('PDF Number (-)'); title(sprintf('Prince coalescence, Wang breakage'));
        grid on; grid minor;    axis square;
        legend('PBM', 'Experiments');
        xlim([0,10])

    end
        

    %
    x = 1;



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