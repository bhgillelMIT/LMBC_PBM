function BreakageTest(params)

%% Setup

    %Settings
    detailed = false;

    %Print update
    fprintf('-- Running Breakage Tests.\n');

    %Load colors
    colors = PlotColors();

    lw = 2;
    fs = 18;


%% Test Luo Svendsen (1997)

if strcmp(params.break.model, 'Luo_Svendson_1996')
    %Define test cases 
    eps = [0.5, 1];
    ds = [0.003, 0.006];

    %Specify line specs
    line_specs = {'k-', 'k--', 'b-', 'b--', 'm-', 'm--'};

    %Create figure
    luofig = figure('Name', 'Luo Svendsen Breakage Test', 'Units', 'normalized', 'Outerposition', [0 0 1 1]);
    subplot(1,2,1);

    %Iterate throguh test cases
    it = 1;
    for ie = 1:length(eps)
        for id = 1:length(ds)

            %Pull test values
            turb_eps = eps(ie);
            d = ds(id);

            %Calculate kolmogorov length scale
            lambda_komogorov = ((params.nus.^3)./(turb_eps+1E-8)).^0.25; %m
            lambda_min = 31.4 * lambda_komogorov; %Minimum eddy diameter to break a bubble - integrate up to bubble diameter 
        
            %Update params
            paramsin = params;
            paramsin.turb.eps = turb_eps;

            %Calculate quantities
            
            [b_eddy, beta_eddy] = BreakageLuoSvendson(1, 1, d, 1, lambda_min, params.break.N_lambdas, paramsin);  


               
            %  %[beta_eddy, beta_ratio] = BreakageInterpolate(d, eps, iz, params);
            % 
            % 
            % 
            % %Calculate breakage rate based on Luo and Svendson 1996
            %  [b_eddy, beta, int_ratio] = BreakageLuoSvendson(1, 10, d, 1, lambda_min, params.break.N_lambdas, params);

            %Plot BSDs
            plot(params.fvs_norm_all, beta_eddy, line_specs{it}, 'LineWidth', 1.5); hold on;

            it = it + 1;
            

        end
    end

    %Update plot aesthetics
    grid on; grid minor; axis square;
    xlabel('Breakage Fraction'); ylabel('Probability');
    title('Luo and Svendson Breakage Test');
    legend('eps = 0.5 m^2/s^3, d = 3 mm', 'eps = 0.5 m^2/s^3, d = 6 mm', 'eps = 1 m^2/s^3, d = 3 mm', 'eps = 1 m^2/s^3, d = 6 mm', 'location', 'northeast');
    set(gca, 'FontSize', fs, 'FontWeight', 'bold');
    ylim([0, 3])
    x = 1;
end
%% Test Wang et al. (2003)

params.break.model = 'Wang_2005';

if strcmp(params.break.model, 'Wang_2005') & false
    epss = [0.25, 0.5, 1.0, 2.0];
    ds = 0.001:0.001:0.01;
    alpha = 0.05;

    %Fill params.alpha_g
    Ns_cell = ones(size(params.cellinds));
    params.alpha_g = alpha .* ones(params.Nz, 1);
    iz = 1;
    
    bd_act_norm = zeros(length(epss), length(ds));
    bd_norm = zeros(length(epss), length(ds));
    errs = bd_norm;

    for ie = 1:length(epss)
        for id = 1:length(ds)

            %Pull test values
            eps = epss(ie);
            d = ds(id);
        
            %Calculate quantities
             [beta_eddy, beta_ratio] = BreakageInterpolate(d, eps, iz, params);

            
            %beta_ratio = params.break.beta_ratio{1}(d, eps);
                %b_eddy = params.break.b_eddy{iz}(d, params.turb.eps(iz));
            %beta_eddy = params.break.beta{1}(d, eps, params.fvs_norm_all);
            %beta_eddy = beta_eddy(:);
            if any(isnan(beta_eddy))
                beta_eddy = zeros(size(beta_eddy));
                beta_ratio = 0;
            end

            %Calculate test quantity
            params.turb.eps = eps;
            N = 1 .* ones(1, params.Nms);
            b_eddy = BreakageEddySimple(1, 10, d, N, beta_ratio, beta_eddy, params);
            bd_norm(ie, id) = b_eddy./((1-alpha)*N(1));


            %Calculate kolmogorov length scale
            lambda_komogorov = ((params.nus.^3)./(eps+1E-8)).^0.25; %m
            lambda_min = 31.4 * lambda_komogorov; %Minimum eddy diameter to break a bubble - integrate up to bubble diameter 


            %Calculate quantity more rigorously
            if detailed
                b_eddy_act = BreakageEddyAlt(1, 10, d, N, lambda_min, params.break.N_lambdas, params);
                err = abs(b_eddy - b_eddy_act)./b_eddy_act;
                if err > 0.0001
                    ratio = b_eddy_act/b_eddy;
                    x = 1;
                end
                errs(ie, id) = err;

                bd_act_norm(ie, id) = b_eddy_act./((1-alpha)*N(1));

            else
                errs(ie, id) = 0;
            end

        end
    end



    %Plot results
    %figure('Name', 'Breakage Tests', 'Units', 'normalized', 'Outerposition', [0 0 1 1]);

    subplot(1,2,2);
    plot(ds, bd_norm(1,:), 'k-','LineWidth', lw, 'Color', colors.trueblue); hold on;
    plot(ds, bd_norm(2,:), 'k-','LineWidth', lw, 'Color', colors.hydrogen); hold on;
    plot(ds, bd_norm(3,:), 'k-','LineWidth', lw, 'Color', colors.truered); hold on;
    plot(ds, bd_norm(4,:), 'k-','LineWidth', lw, 'Color', colors.ChiliRed); hold on;

    plot(ds, bd_act_norm(1,:), 'k--','LineWidth', lw, 'Color', colors.trueblue); hold on;
    plot(ds, bd_act_norm(2,:), 'k--','LineWidth', lw, 'Color', colors.hydrogen); hold on;
    plot(ds, bd_act_norm(3,:), 'k--','LineWidth', lw, 'Color', colors.truered); hold on;
    plot(ds, bd_act_norm(4,:), 'k--','LineWidth', lw, 'Color', colors.ChiliRed); hold on;

    grid on; grid minor; axis square;
    xlabel('Diameter (m)'); ylabel('$\frac{b(d)}{[(1-\alpha_d)N]}$', 'Interpreter', 'latex' );
    title('Normalized Breakage Rate');

    %Create legend
    h = zeros(1,6);
    h(1) = plot(nan, nan, 'k-', 'LineWidth', lw); hold on;
    h(2) = plot(nan, nan, 'k-', 'LineWidth', lw); hold on;
    h(3) = plot(nan, nan, 'k-', 'LineWidth', lw, 'Color', colors.trueblue); hold on;
    h(4) = plot(nan, nan, 'k-', 'LineWidth', lw, 'Color', colors.hydrogen); hold on;
    h(5) = plot(nan, nan, 'k-', 'LineWidth', lw, 'Color', colors.truered); hold on;
    h(6) = plot(nan, nan, 'k-', 'LineWidth', lw, 'Color', colors.ChiliRed); hold on;
    legend(h, 'Interpolated', 'Actual', '\epsilon = 0.25 m^2/s^3', '\epsilon = 0.5 m^2/s^3', '\epsilon = 1.0 m^2/s^3', '\epsilon = 2.0 m^2/s^3', 'location', 'northwest');
    %legend('\epsilon = 0.25 m^2/s^3', '\epsilon = 0.5 m^2/s^3', '\epsilon = 1.0 m^2/s^3', '\epsilon = 2.0 m^2/s^3', 'location', 'northwest');
    set(gca, 'FontSize', fs, 'FontWeight', 'bold');
end
    
%% BSD Tests

    %Test cases
    subplot(1,2,2);
    %figure();
    eps_ref = 1;
    ds_ref = [0.0015, 0.002, 0.003, 0.006];
    BSDs = zeros(length(ds_ref), length(params.fvs_norm_all));
    iz = 1;
    for it_debug = 1:length(ds_ref)

        %Calculated interpolated value
        d = ds_ref(it_debug);
        [beta_eddy, beta_ratio] = BreakageInterpolate(d, eps_ref, iz, params); %params.break.funcs{iz}.betas{iz}(ds_ref(it_debug), eps_ref, params.fvs_norm_all);
        BSDs(it_debug,:) = beta_eddy;
        %Renormalize interpolated value

        %Plot
        %subplot(1,2,1);
        plot(params.fvs_norm_all, BSDs(it_debug,:), 'LineWidth', 1.5); hold on;

        % %Calculate manually
        % turb_eps = eps_ref;
        % params.turb.eps(:) = turb_eps;
        % lambda_komogorov = ((params.nus.^3)./(turb_eps+1E-8)).^0.25; %m
        % lambda_min = 31.4 * lambda_komogorov; %Minimum eddy diameter to break a bubble - integrate up to bubble diameter 
        % [b_eddy, beta, int_ratio] = BreakageEddyAlt(iz, id, d, Ns_cell, lambda_min, params.break.N_lambdas, params);
        % subplot(1,2,2);
        % plot(params.fvs_norm_all, beta, 'LineWidth', 1.5); hold on;

    end


    %Plot aesthetics
    grid on; grid minor; axis square;
    xlabel('Breakage Fraction'); ylabel('Probability');
    title('Bubble Size Distributions');
    legend('d = 1.5 mm', 'd = 2 mm', 'd = 3 mm', 'd = 6 mm', 'location', 'north');
    set(gca, 'FontSize', fs, 'FontWeight', 'bold');


    %Review and close
    %close;

    %Pause to render figure
    pause(1);


%% Compare Wang and Luo models for the same conditions 

    %Create new figure
    figure('Name', 'Wang vs Luo Breakage', 'Units', 'normalized', 'Outerposition', [0 0 1 1]);

    %Define inputs
    epss = [0.5, 2];
    ds = [0.0025, 0.01];
    alpha = 0.05;

    %Define line specs
    linespecs = {'-', ':', '-.', '--'};

    %Iterate through the test cases 
     for ie = 1:length(epss)

        subplot(1,length(epss), ie);
        for id = 1:length(ds)

            %Pull test values
            turb_eps = epss(ie);
            d = ds(id);

            %Calculate kolmogorov length scale
            lambda_komogorov = ((params.nus.^3)./(turb_eps+1E-8)).^0.25; %m
            lambda_min = 31.4 * lambda_komogorov; %Minimum eddy diameter to break a bubble - integrate up to bubble diameter 
        
            %Update params
            paramsin = params;
            paramsin.break.model = 'Wang_2005';
            paramsin.turb.eps = turb_eps;


            %Calculate Wang result
            [beta_wang, ~, ~] = BreakageInterpolate(d, turb_eps, 1, paramsin);


            %Calculate Luo result
            [~, beta_luo] = BreakageLuoSvendson(1, 1, d, 1, lambda_min, params.break.N_lambdas, paramsin);
            
            %Specify line specs
            wang_spec = ['k', linespecs{id}];
            luo_spec = ['b', linespecs{id}];

            %Plot the results
            plot(params.fvs_norm_all, beta_wang, wang_spec, 'LineWidth', 1.5); hold on;
            plot(params.fvs_norm_all, beta_luo, luo_spec, 'LineWidth', 1.5); 





        end

        %Define plot aesthetics
        grid on; grid minor; axis square;
        xlabel('Breakage Fraction'); ylabel('Probability');
        title('Wang vs Luo Breakage (\epsilon = ' + string(epss(ie)) + ' m^2/s^3)');
        ylim([0,20])
        set(gca, 'FontSize', fs, 'FontWeight', 'bold');

        %Create legend
        h = zeros(1,4);
        h(1) = plot(nan, nan, 'k.', 'MarkerSize', 18); hold on;
        h(2) = plot(nan, nan, 'b.', 'MarkerSize', 18); hold on;
        h(3) = plot(nan, nan, 'k-', 'LineWidth', 1.5, 'Color', [0.6, 0.6, 0.6]); hold on;
        h(4) = plot(nan, nan, 'k:', 'LineWidth', 1.5, 'Color', [0.6, 0.6, 0.6]); hold on;
        %h(5) = plot(nan, nan, 'k-.', 'LineWidth', 1.5, 'Color', [0.6, 0.6, 0.6]); hold on;
        legend(h, 'Wang', 'Luo', 'd = 2.5 mm', 'd = 10 mm', 'location', 'northeast');


        
    end

   

    

    
    


end


%Calculate the breakage rate from interpolated values of 
function b_eddy = BreakageEddySimple(iz, im, d, Ns_cell, beta_ratio, beta_eddy, params)

    %Debug cases
    if params.break.debug
        %d = 0.003;
        t_start.total = cputime;
        fprintf('-Eddy Breakage: iz = %i; im = %i; d = %0.5f m; m_i = %0.3e;\n', iz, im, d, params.mms(im));
        x = 1;
    end

    %Track iterations
    it_total = 1;

    %Calculate integral
    b_fvd = params.break.bfd_zero;
    int_lambda_fvs = beta_ratio .* beta_eddy(1:length(params.fvs_norm));
    b_fvd = 0.923 .* (1 - params.alpha_g(iz)) .* Ns_cell(im) .* params.turb.eps(iz).^(1/3)  .* int_lambda_fvs;

    % for iv = 2:length(params.fvs_norm)
    % 
    %     %Calulcate 
    %     int_lambda_fv = beta_ratio .* beta_eddy(iv);
    %     b_fvd(iv) = 0.923 .* (1 - params.alpha_g(iz)) .* Ns_cell(im) .* params.turb.eps(iz).^(1/3) .* int_lambda_fv;
    % 
    %     %Update iteration counter
    %     it_total = it_total + 1;
    % 
    % end

    %Calculate overall rate
    switch params.break.int_method
        case 'gauss'
            int_func = @(x) interp1(params.fvs_norm, b_fvd, x);
            b_eddy = quadgk(int_func, 0, 1);
        case 'trapz'
    	    b_eddy = trapz(params.fvs_norm, b_fvd);
        otherwise
            warning('Invalid Integration Method. Defaulting to trapz.');
            b_eddy = trapz(params.fvs_norm, b_fvd);
    end

    %Debug plot
    if params.break.debug
        t_end = cputime;
        t_req = t_end - t_start.total;
        % figure();
        % plot(params.fvs_norm_all, beta_eddy);
        % close
        fprintf('-- Eddy Breakage time = %0.8f; Total its = %d; \n', t_req, it_total)
   end


    

end