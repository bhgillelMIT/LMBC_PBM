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



end