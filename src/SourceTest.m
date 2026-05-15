function SourceTest(params)

    %Settings
    debug = true;
    lw = 2;
    fs = 16;
    disttype = 'number'; %'number' or 'mass' - determines whether the distribution of source terms is normalized to number or mass

    %Load colors


    %Create debug figure
    if debug
        debugfig = figure();
    end

    %Make params global
    global params
    params.src.its  = 1;

    %Copy params
    paramsin = params;
    paramsin.src.solve_eps = false;
    paramsin.src.solve_alphag = false;
    paramsin.src.alphag_manual = 0.1;
    paramsin.sol.single_layer = true;


    %Update input values to simplify analysis
    Nz_actual = params.Nz;
    src_its_actual = 0;
    paramsin.Nz = 1;
    paramsin.src.its = 1;
    paramsin.turb.k = 0.1;

    %Define parameters for tests
    epss = [0.2, 1, 5]; %m^2/s^3 - Turbulent energy dissipation rates to test
    alpha_g = 0.1;

    %Calculate initial distribution based on the gas holdup, assuming a uniform distribution
    switch disttype
        case 'number'
            yin = 30000*ones(1, params.Nms); %yin = params.N_dot_os./trapz(params.N_dot_os);
        case 'mass'
            Vgas = alpha_g;
            Vpdf = ones(1, params.Nms)./(params.mms(end)); %pdf normalized to mass
            Vgass = Vpdf./params.Nms .* ones(1, params.Nms);
            Ns = Vgass./params.Vms; %Numeric density 
            yin = Ns;
    end

    

    %Storage matricies
    cadds = zeros(length(epss), params.Nms); %1/s - Total coalescence rate - additive term
    csubs = zeros(length(epss), params.Nms); %1/s - Total coalescence rate - subtractive term
    cts = zeros(length(epss), params.Nms); %1/s - Coalescence rate due to turbulent eddies
    cws = zeros(length(epss), params.Nms); %1/s - Coalescence rate due to wake entrainment
    crs = zeros(length(epss), params.Nms); %1/s - Coalescence rate due to rise velocity differences
    badds = zeros(length(epss), params.Nms); %1/s - Total breakage rate - additive term
    bsubs = zeros(length(epss), params.Nms); %1/s - Total breakage rate - subtractive term
    bts = zeros(length(epss), params.Nms); %1/s - Breakage rate due to turbulent eddies
    bss = zeros(length(epss), params.Nms); %1/s - Breakage rate due to surface instability
    hs = zeros(length(epss), params.Nms); %1/s - Net source term (coalescence - breakage)


    %Iterate through each value of eps and calculate the source term for each size/spatial cell
    for ie = 1:length(epss)

        %Pull value of eps and d
        paramsin.src.eps_manual = epss(ie);
        paramsin.turb.eps = epss(ie);


        %Calculate local properties based on this value of eps
        paramsin = CalcLocalProperties(yin, paramsin, true);

        %Test Wang 2005 models
        paramsin.break.model = 'Wang_2005';
        paramsin.coalesce.model = 'Wang_2005';



        %Calculate coalescence - turbulence only - yin = 1000 
        paramsin.coalesce.eddy = true;
        paramsin.coalesce.wake = false;
        paramsin.coalesce.rise = false;
        [cadd_eddy, csub_eddy, ~] = Coalescence(yin, paramsin, true);

        % %Calculate coalescence - turbulence only - yin = 10000
        % paramsin.coalesce.eddy = true;
        % paramsin.coalesce.wake = false;
        % paramsin.coalesce.rise = false;
        % [cadd_eddy_10000, csub_eddy, ~] = Coalescence(10*yin, paramsin, true);
        % 
        % %Calculate coalescence - turbulence only - yin = 25000
        % paramsin.coalesce.eddy = true;
        % paramsin.coalesce.wake = false;
        % paramsin.coalesce.rise = false;
        % [cadd_eddy_25000, csub_eddy, ~] = Coalescence(25*yin, paramsin, true);

        %Calculate coalescence - wake only - yin = 1000
        paramsin.coalesce.eddy = false;
        paramsin.coalesce.wake = true;
        paramsin.coalesce.rise = false;
        [cadd_wake, csub_wake, ~] = Coalescence(yin, paramsin, true);

        %Calculate coalescence - rise velocity only - yin = 1000
        paramsin.coalesce.eddy = false;
        paramsin.coalesce.wake = false;
        paramsin.coalesce.rise = true;
        [cadd_rise, csub_rise, ~] = Coalescence(yin, paramsin, true);

        %Calculate coalescence - all - yin = 1000
        paramsin.coalesce.eddy = true;
        paramsin.coalesce.wake = true;
        paramsin.coalesce.rise = true;
        [cadd_wang, csub_wang, ~] = Coalescence(yin, paramsin, true);

        %Calculate coalescence - all - yin = 1000
        paramsin.coalesce.eddy = true;
        paramsin.coalesce.wake = false;
        paramsin.coalesce.rise = true;
        [cadd_wang_turb_rise, csub_wang_turb_rise, ~] = Coalescence(yin, paramsin, true);

        %Calculate coalescence - prince and blanch - yin = 1000
        paramsin.coalesce.model = 'PrinceBlanch_1990';
        paramsin.coalesce.eddy = true;
        [cadd_pb, csub_pb, ~] = Coalescence(yin, paramsin, true);

        %Calculate breakage - turbulence only - yin = 1000
        paramsin.break.eddy = true;
        paramsin.break.surface = true;
        [badd_wang, bsub_wang, bmats] = Breakage(yin, paramsin, true);

        %Calculate breakage - luo and svendson - yin = 1000
        paramsin.break.model = 'Luo_Svendson_1996';
        paramsin.break.eddy = true;
        [badd_luo, bsub_luo, bmats] = Breakage(yin, paramsin, true);

        %Calculate net source terms
        h_wang = cadd_wang(:) - csub_wang(:) + badd_wang(:) - bsub_wang(:);
        h_prince_luo = cadd_pb(:) - csub_pb(:) + badd_luo(:) - bsub_luo(:);

        %Pause and plot results for this eps
        if debug

            %Select figure
            figure()

            %Plot coalescence rates
            subplot(1,3,1);
            semilogy(params.dms, cadd_eddy, 'b-', 'LineWidth', lw); hold on;
            semilogy(params.dms, csub_eddy, 'b--', 'LineWidth', lw);
            semilogy(params.dms, cadd_wake, 'r-', 'LineWidth', lw);
            semilogy(params.dms, csub_wake, 'r--', 'LineWidth', lw);
            semilogy(params.dms, cadd_rise, 'm-', 'LineWidth', lw);
            semilogy(params.dms, csub_rise, 'm--', 'LineWidth', lw);
            semilogy(params.dms, cadd_pb, 'c-', 'LineWidth', lw);
            semilogy(params.dms, csub_pb, 'c--', 'LineWidth', lw);
            
            xlabel('Diameter (m)'); ylabel('Coalescence Rate (1/s)'); title(['Coalescence (\epsilon = ', num2str(epss(ie)), ')']);
            grid on; grid minor; axis square;  
            legend('cadd - eddy', 'csub - eddy', 'cadd - wake', 'csub - wake', 'cadd - rise', 'csub - rise', 'cadd - pb', 'csub - pb', 'location', 'southeast');
            
            ylim([1E-3, 1E8])

            subplot(1,3,2);
            semilogy(params.dms, badd_wang, 'b-', 'LineWidth', lw); hold on;
            semilogy(params.dms, bsub_wang, 'b--', 'LineWidth', lw);
            semilogy(params.dms, badd_luo, 'r-', 'LineWidth', lw);
            semilogy(params.dms, bsub_luo, 'r--', 'LineWidth', lw);
            xlabel('Diameter (m)'); ylabel('Breakage Rate (1/s)'); title(['Breakage (\epsilon = ', num2str(epss(ie)), ')']);
            grid on; grid minor; axis square;
            legend('badd - wang', 'bsub - wang', 'badd - luo', 'bsub - luo', 'location', 'southeast');
            ylim([1E-3, 1E8])

            %Plot total source term
            subplot(1,3,3);
            plot(params.dms, h_wang, 'b-', 'LineWidth', lw); hold on;
            plot(params.dms, h_prince_luo, 'r-', 'LineWidth', lw);
            xlabel('Diameter (m)'); ylabel('Net Source Term (1/s)'); 
            title(['Net Source Term (\epsilon = ', num2str(epss(ie)), ')']);
            grid on; grid minor; axis square;
            legend('Wang 2005', 'Prince/Luo', 'location', 'southeast');
            %ylim([1E-3, 1E8])
            
            

        end

        

        %Plot a comparison of wang and luo models
        figure();

        semilogy(params.dms, cadd_eddy, 'b-', 'LineWidth', lw, 'Color', [20/255, 52/255, 164/255]); hold on;
        semilogy(params.dms, cadd_wang_turb_rise, 'b-', 'LineWidth', lw, 'Color', [0, 150/255, 1]); hold on;
        semilogy(params.dms, cadd_wang, 'b-', 'LineWidth', lw, 'Color', [0,1,1])
        semilogy(params.dms, cadd_pb, 'c-', 'LineWidth', lw, 'Color', [0.8, 0.2, 0.1]); hold on;
        semilogy(params.dms, csub_eddy, 'b--', 'LineWidth', lw, 'Color', [20/255, 52/255, 164/255]); hold on;
        semilogy(params.dms, csub_wang_turb_rise, 'r--', 'LineWidth', lw, 'Color', [0, 150/255, 1]); 
        semilogy(params.dms, csub_wang, 'r--', 'LineWidth', lw, 'Color', [0,1,1]);
        semilogy(params.dms, csub_pb, 'c--', 'LineWidth', lw, 'Color', [1, 0.2, 0.1]); hold on;
        xlabel('Diameter (m)'); ylabel('Coalescence Rate (1/s)'); 
        title(['Coalescence Comparison (\epsilon = ', num2str(epss(ie)), ')']);
        grid on; grid minor; axis square;

    

        legend('c^+_{wang} - turb', 'c^+_{wang} - turb + rise', 'c^+_{wang} - turb + rise + wake',  'c^+_{prince} - turb + rise',...
            'c^-_{wang} - turb', 'c^-_{wang} - turb + rise','c^-_{wang} - turb + rise + wake', 'c^-_{prince} - turb + rise',...
            'location', 'southeast');
        ylim([1E-3, 1E8]) 
        set(gca, 'FontSize', fs);


        %Log results



        x = 1;




        % %Log results
        % cadds(ie,:) = cadd';
        % badds(ie,:) = badd';
        % hs(ie,:) = h';



        

    end


%% Wang 2005 - Comparison of Coalescence rates at different 

    %Settings
    disttype = 'number';

    %Create new figure
    wangfig = figure();

    %Test conditions
    epss = [0.2, 1, 5];
    alpha_gs = [0.01, 0.1, 0.3];

    %Linespecs
    colors = [0.1, 0.3, 1;
              0.5, 0.3, 0.5;
              1, 0.3, 0.1];
    linespecs = {'-', '--', '-.'};

    %Iterate through gas holdups
    for ia = 1:length(alpha_gs)

        %Pull gas holdup
        alpha_g = alpha_gs(ia);

        %Calculate initial distribution based on the gas holdup, assuming a uniform distribution
        switch disttype
            case 'number'
                Vgas = alpha_g;
                Vgas_func = @(N) sum(N .* ones(1, params.Nms) .* params.Vms) - Vgas;
                N = fzero(Vgas_func, 10);
                
                yin = N*ones(1, params.Nms); %yin = params.N_dot_os./trapz(params.N_dot_os);
            case 'mass'
                Vgas = alpha_g;
                Vpdf = ones(1, params.Nms)./(params.mms(end)); %pdf normalized to mass
                Vgass = Vpdf./params.Nms .* ones(1, params.Nms);
                Ns = Vgass./params.Vms; %Numeric density 
                yin = Ns;
        end

        %Iterate through eps values
        for ie = 1:length(epss)

            %Define linespec
            linespec = ['k', linespecs{ia}];
            linecolor = colors(ie, :);

            %Update params
            paramsin.src.eps_manual = epss(ie);
            paramsin.turb.eps = epss(ie);
            paramsin.src.alphag_manual = alpha_g;
            paramsin.alpha_g = alpha_g;

            %Calculate local properties based on this value of eps and alpha_g
            paramsin = CalcLocalProperties(yin, paramsin, true);

            %Calculate coalescence - turbulence only
            paramsin.coalesce.eddy = true;
            paramsin.coalesce.wake = false;
            paramsin.coalesce.rise = false;
            [cadd_eddy, csub_eddy, ~] = Coalescence(yin, paramsin, true);

            %Calculate coalescence - wake only
            paramsin.coalesce.eddy = false;
            paramsin.coalesce.wake = true;
            paramsin.coalesce.rise = false;
            [cadd_wake, csub_wake, ~] = Coalescence(yin, paramsin, true);

            %Calculate coalescence - rise velocity only
            paramsin.coalesce.eddy = false;
            paramsin.coalesce.wake = false;
            paramsin.coalesce.rise = true;
            [cadd_rise, csub_rise, ~] = Coalescence(yin, paramsin, true);

            %Calculate caolescence - all
            paramsin.coalesce.eddy = true;
            paramsin.coalesce.wake = true;
            paramsin.coalesce.rise = true;
            [cadd_all, csub_all, ~] = Coalescence(yin, paramsin, true);

            %Plot result
            semilogy(params.dms, cadd_all./N, linespec, 'LineWidth', lw, 'Color', linecolor); hold on;

            


        end
    end

    %Define plot aesthetics
    xlabel('Diameter (m)'); ylabel('Coalescence Rate (1/s)'); 
    title('Coalescence Rate from Wang 2005');
    grid on; grid minor; axis square;
    ylim([1E-3, 1E2]);

    %Define legend
    h = zeros(1, 6);
    h(1) = plot(nan, nan, 'k-', 'LineWidth', lw); hold on;
    h(2) = plot(nan, nan, 'k--', 'LineWidth', lw); hold on;
    h(3) = plot(nan, nan, 'k-.', 'LineWidth', lw); hold on;
    h(4) = plot(nan, nan, 'k-', 'LineWidth', lw, 'Color', colors(1,:)); hold on;
    h(5) = plot(nan, nan, 'k-', 'LineWidth', lw, 'Color', colors(2,:)); hold on;
    h(6) = plot(nan, nan, 'k-', 'LineWidth', lw, 'Color', colors(3,:)); hold on;
    legend(h, '\alpha_g = 0.01', '\alpha_g = 0.1', '\alpha_g = 0.3', '\epsilon = 0.2 m^2/s^3', '\epsilon = 1 m^2/s^3', '\epsilon = 5 m^2/s^3', 'location', 'southeast');

    %Update font size
    set(gca, 'FontSize', fs);

    %Itera


    %Restore values for params
    params.src.its  = 0;




end