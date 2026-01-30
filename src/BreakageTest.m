function BreakageTest(params)

%% Setup

    %Print update
    fprintf('-- Running Breakage Tests.\n');

    %Load colors
    colors = PlotColors();

    lw = 2;
    fs = 18;


%% Test Wang et al. (2003)

    epss = [0.25, 0.5, 1.0, 2.0];
    ds = 0.001:0.001:0.01;
    alpha = 0.05;

    %Fill params.alpha_g
    Ns_cell = ones(size(params.cellinds));
    params.alpha_g = alpha .* ones(params.Nz, 1);
    iz = 1;
    

    bd_norm = zeros(length(epss), length(ds));

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

        end
    end



    %Plot results
    %figure('Name', 'Breakage Tests', 'Units', 'normalized', 'Outerposition', [0 0 1 1]);

    subplot(1,2,1);
    plot(ds, bd_norm(1,:), 'k-','LineWidth', lw, 'Color', colors.trueblue); hold on;
    plot(ds, bd_norm(2,:), 'k-','LineWidth', lw, 'Color', colors.hydrogen); hold on;
    plot(ds, bd_norm(3,:), 'k-','LineWidth', lw, 'Color', colors.truered); hold on;
    plot(ds, bd_norm(4,:), 'k-','LineWidth', lw, 'Color', colors.ChiliRed); hold on;

    grid on; grid minor; axis square;
    xlabel('Diameter (m)'); ylabel('$\frac{b(d)}{[(1-\alpha_d)N]}$', 'Interpreter', 'latex' );
    title('Normalized Breakage Rate');
    legend('\epsilon = 0.25 m^2/s^3', '\epsilon = 0.5 m^2/s^3', '\epsilon = 1.0 m^2/s^3', '\epsilon = 2.0 m^2/s^3', 'location', 'northwest');
    set(gca, 'FontSize', fs, 'FontWeight', 'bold');

    
%% BSD Tests

    %Test cases
    %subplot(1,2,2);
    figure();
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
        subplot(1,2,1);
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