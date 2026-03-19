function params = CalculateBSDs(params)

    %Calculate range of gas velocities
    u_spf_min = params.u_spf_orifice;
    u_spf_max = params.u_spf_surface;
    u_spfs = linspace(0.9*u_spf_min, 1.1*u_spf_max, params.break.N_u_spfs);

    %Calculate range of bubble diameters
    V_min = (params.nms(1) .* params.R .* params.T_mu_i)./params.p_orifice;
    V_max = (2.*params.nms(end) .* params.R .* params.T_liq)./params.p_surf;
    d_min = (6.*V_min/pi).^(1/3); 
    d_max = (6.*V_max/pi).^(1/3);

    ds = logspace(log10(params.break.d_range(1)), log10(params.break.d_range(2)), params.break.N_ds);

    %Calculate range of turbulent kinetic energy dissipation (epsilon)
    turb_epss = [0, logspace(log10(0.1), log10(params.break.eps_range(2)), params.break.N_u_spfs-1)]; %params.g * u_spfs';

    
    %Determine boundaries for interpolation functions
    N_funcs = params.break.N_d_interp * params.break.N_eps_interp;
    if length(turb_epss) > params.break.N_eps_interp
        epss_inds = 1:length(turb_epss)/params.break.N_eps_interp:length(turb_epss);
        epss_inds(end+1) = length(turb_epss);
        epss_inds = round(epss_inds);
    else
        error('Fewer epsilon points than subdivisions.')
    end

    if length(ds) > params.break.N_d_interp
        ds_inds = 1:length(ds)/params.break.N_d_interp:length(ds);
        ds_inds(end+1) = length(ds);
        ds_inds = round(ds_inds);
    else
        error('Fewer d points than subdivisions.');
    end
    ds_total = ds;
    epss_total = turb_epss;
    

    




    %Determine number of spatial cells to iterate through
    if params.isothermal
        Nz = 1;
    else
        Nz = params.Nz;
    end

    %Allocate storage arrays for total dataset
    ds_mg_total = zeros(params.break.N_ds, params.break.N_eps, length(params.fvs_norm_all));
    ds_mg2_total = zeros(params.break.N_ds, params.break.N_eps);
    turb_epss_mg_total = zeros(params.break.N_ds, params.break.N_eps, length(params.fvs_norm_all));
    turb_epss_mg2_total = zeros(params.break.N_ds, params.break.N_eps);
    fvs_mg_total = zeros(params.break.N_ds, params.break.N_eps, length(params.fvs_norm_all));
    betas_total = zeros(params.break.N_ds, params.break.N_eps, length(params.fvs_norm_all));
    bs_ints_ratio_total = zeros(params.break.N_ds, params.break.N_eps);
    bs_total = bs_ints_ratio_total;


    %Iterate through spatial cells
    %params.bs = cell(1, params.Nz);
    params.break.funcs_d_range = zeros(N_funcs, 2);
    params.break.funcs_eps_range = zeros(N_funcs, 2);
    params.break.funcs = cell(1, Nz);
    params.beta_ratio = cell(1, params.Nz);
    params.betas = cell(1, params.Nz);
    for iz = 1:Nz

        %Allocate storage cell array
        params.break.funcs{iz} = cell(1, N_funcs);

        func_it = 1;
        %Iterate through subdivisions
        for ie = 1:params.break.N_eps_interp
            for id1 = 1:params.break.N_d_interp

                %Pull bounds
                ds = ds_total(ds_inds(id1):ds_inds(id1+1));
                turb_epss = epss_total(epss_inds(ie):epss_inds(ie+1));

                %Create mesh grids
                P = [2, 1, 3];
                [ds_mg, turb_epss_mg, fvs_mg] = meshgrid(ds, turb_epss, params.fvs_norm_all);
                ds_mg = permute(ds_mg, P); turb_epss_mg = permute(turb_epss_mg, P);
                fvs_mg = permute(fvs_mg, P);
                [eps_mg2, ds_mg2] = meshgrid(turb_epss, ds);
                betas = zeros(size(ds_mg));
                bs = zeros(length(ds), length(turb_epss));
                bs_ints_ratio = bs; %zeros(size(ds_mg));
                %bs_ints_ratio = bs_ints_ratio(:,:,1:length(params.fvs_norm));
        
                %Pull numeric densities in this spatial celll
                z = params.cents_y(iz);
                cellinds = ((iz-1).*params.Nms + 1):1:(iz.*params.Nms);
                %Ns = ones(size(cellinds));
                Ns_cell = ones(size(cellinds));
        
                %Pull local fluid properties
                uL = 0;
        
                %Calculate critical breakup diameter
                d_crit = CalcCritDiameter(iz, params); %m - critical diameter for instability breakage
        
                %Create storage matrix



                %VARY RESOLUTION OF fvs norm based on the number of size groups
                %below it to reduce computational cost
        
                %Iterate through size groups
        
                for id = 1:length(ds)
            
                    %Pull bubble size
                    %mi = params.mms(im); %Current represenative mass
                    %V = (params.nms(im) .* (1 + params.X_mu(iz, im)) .* params.R .* params.T_mu(iz,im))./params.p_func(z);
                    d = ds(id); %(6.*V/pi).^(1/3); d = 0.003; 

                    %u = params.uzs(cellinds(im));
        
        
        
                    %Iterate through possible epsilons
                    paramsin = params;
                    N_lambdas = params.break.N_lambdas;
                    for ie2 = 1:length(turb_epss)
            
                        %Pull epsilon
                        turb_eps = turb_epss(ie2);
                        paramsin.turb.eps(:) = turb_eps;
            
                        %Calculate kolmogorov length scale
                        lambda_komogorov = ((paramsin.nus.^3)./(turb_eps+1E-8)).^0.25; %m
                        lambda_min = 31.4 * lambda_komogorov; %Minimum eddy diameter to break a bubble - integrate up to bubble diameter 
        
                        %Calculate minimum lambda for breakage
                        
            
                        %Calculate breakage rate and BSD
                        switch params.break.model
                            case 'Luo_Svendson_1996'
                                [b_eddy, beta, int_ratio] = BreakageLuoSvendson(iz, id, d, Ns_cell, lambda_min, N_lambdas, paramsin);
                            case 'Wang_2005'
                                [b_eddy, beta, int_ratio] = BreakageEddyAlt(iz, id, d, Ns_cell, lambda_min, N_lambdas, paramsin);
                            otherwise
                                error('Breakage model not recognized. Options are "Luo_Svendson_1996", and "Wang_2005".');
                        end

                        %Store results
                        bs_ints_ratio(id, ie2) = int_ratio;
                        bs(id, ie2) = b_eddy; %Breakage rate
                        betas(id, ie2, :) = beta; %Bubble size distribution
                        
            
                    end
        
                    
                end

                %Store data - ind 1 = d; ind 2 = eps; ind3 = fv
                
        

                if id1 == 1 & ie == 1
                    ds_mg_inds = ds_inds(id1):1:ds_inds(id1+1);
                    epss_mg_inds = epss_inds(ie):1:epss_inds(ie+1);        
                else
                    ds_mg_inds = ds_inds(id1):1:ds_inds(id1+1);
                    epss_mg_inds = epss_inds(ie):1:epss_inds(ie+1); 

                    x = 1;
                end

                %Store results
                ds_mg_total(ds_mg_inds, epss_mg_inds, :) = ds_mg;
                ds_mg2_total(ds_mg_inds, epss_mg_inds) = ds_mg2;
                turb_epss_mg_total(ds_mg_inds, epss_mg_inds, :) = turb_epss_mg;
                turb_epss_mg2_total(ds_mg_inds, epss_mg_inds) = eps_mg2;
                fvs_mg_total(ds_mg_inds, epss_mg_inds, :) = fvs_mg;
                betas_total(ds_mg_inds, epss_mg_inds, :) = betas;
                bs_ints_ratio_total(ds_mg_inds, epss_mg_inds) = bs_ints_ratio;
                bs_total(ds_mg_inds, epss_mg_inds) = bs;


                %Create and store interpolation function
                %betas = permute(betas, P);
                %F = griddedInterpolant(ds_mg, turb_epss_mg, fvs_mg, betas);

                params.break.funcs{iz}{func_it}.eps_range = [min(turb_epss), max(turb_epss)];
                params.break.funcs{iz}{func_it}.d_range = [min(ds), max(ds)];
                params.break.funcs{iz}{func_it}.b_eddy = @(d, eps) interp2(ds_mg2', eps_mg2', bs, d, eps);
                params.break.funcs{iz}{func_it}.betas = @(d, eps, fv) interpn(ds_mg, turb_epss_mg, fvs_mg, betas, d, eps, fv); %function of fv, d, eps
                params.break.funcs{iz}{func_it}.beta_ratio = @(d, eps) interp2(ds_mg2', eps_mg2', bs_ints_ratio', d, eps); %params.bs{iz} = @(d, eps) interp2(ds_mg2, eps_mg2, bs, d, eps);
                params.break.funcs_eps_range(func_it,:) = [min(turb_epss), max(turb_epss)];
                params.break.funcs_d_range(func_it,:) = [min(ds), max(ds)];


                %
                func_it = func_it + 1;

            end

            



        end

        %Create full cost interpolation functions
        %params.break.func.eps_range = [min(turb_epss), max(turb_epss)];
        %params.break.func.d_range = [min(ds), max(ds)];
        params.break.b_eddy = @(d, eps) interp2(ds_mg2_total', eps_mg2_total', bs_total', d, eps);
        params.break.betas = @(d, eps, fv) interpn(ds_mg_total, turb_epss_mg_total, fvs_mg_total, betas_total, d, eps, fv); %function of fv, d, eps
        params.break.beta_ratio = @(d, eps) interp2(ds_mg2_total', eps_mg2_total', bs_ints_ratio_total', d, eps); %params.bs{iz} = @(d, eps) interp2(ds_mg2, eps_mg2, bs, d, eps);
        
        
        %params.break.funcs_eps_range(func_it,:) = [min(turb_epss), max(turb_epss)];
        %params.break.funcs_d_range(func_it,:) = [min(ds), max(ds)];
            



    end


    % %Test the 
    % if params.debug
    % 
    %     %Create figure
    %     figure();
    % 
    %     %Reference cases
    % 
    %     %Test cases
    %     eps_ref = 1;
    %     ds_ref = [0.0015, 0.002, 0.003, 0.006];
    %     BSDs = zeros(length(ds_ref), length(params.fvs_norm_all));
    %     for it_debug = 1:length(ds_ref)
    % 
    %         %Calculated interpolated value
    %         BSDs(it_debug,:) = params.break.funcs{iz}.betas{iz}(ds_ref(it_debug), eps_ref, params.fvs_norm_all);
    % 
    %         %Renormalize interpolated value
    % 
    % 
    %         %Plot
    %         plot(params.fvs_norm_all, BSDs(it_debug,:), 'LineWidth', 1.5); hold on;
    % 
    %     end
    % 
    % 
    %     %Plot aesthetics
    %     grid on; grid minor; axis square;
    % 
    % 
    %     %Review and close
    %     %close;
    % 
    % 
    % 
    % 
    % end

    %If isothermal, convert to standard format for each z-level
    if params.isothermal
        for i = 2:params.Nz
            for func_it = 1:N_funcs
                % params.betas{i} = params.betas{1};
                %params.bs{2:end} = params.bs{1};
    
                params.break.funcs{i}{func_it}.beta_ratio = params.break.funcs{1}{func_it}.beta_ratio;
                params.break.funcs{i}{func_it}.betas = params.break.funcs{1}{func_it}.betas;
                params.break.funcs{i}{func_it}.eps_range = params.break.funcs{1}{func_it}.eps_range;
                params.break.funcs{i}{func_it}.d_range = params.break.funcs{1}{func_it}.d_range;

            end
        end
    end

    %Consider BSDs as a function of d and epsilon for each spatial cell,
    %which is assumed to have fixed conditions 

    %Save interpolation function
    break_out.funcs = params.break.funcs;
    break_out.d_range =  params.break.funcs_d_range;
    break_out.eps_range = params.break.funcs_eps_range;
    
    %Save file
    switch params.break.model
        case 'Wang_2005'
            outname = sprintf('break_%s_Wang_Nd-%d_Nz-%d_Ne-%d_TL-%d.mat', params.liquid.name, params.break.N_d_interp, params.Nz, params.break.N_eps_interp, round(params.T_liq));
        case 'Luo_Svendson_1996'
            outname = sprintf('break_%s_LuoSvendson_Nd-%d_Nz-%d_Ne-%d_TL-%d.mat', params.liquid.name, params.break.N_d_interp, params.Nz, params.break.N_eps_interp, round(params.T_liq));
        otherwise
            error('Invalid breakage model. Options are "Wang_2005", and "Luo_Svendson_1996".');
    end
            
    outname = sprintf('break_%s_Nd-%d_Nz-%d_Ne-%d_TL-%d.mat', params.liquid.name, params.break.N_ds, params.Nz, params.break.N_u_spfs, round(params.T_liq));
    outpath = [params.folders.break, outname];
    save(outname, 'break_out');

    %Define outputs - to continue 
    params.break.funcs = break_out.funcs;
    params.break.funcs_d_range = break_out.d_range;
    params.break.funcs_eps_range = break_out.eps_range;

    %Check that the functions work
    if true
            
            %Create figure
            figure();
            
            %Test each of the windows
            for i = 1:N_funcs

                %Specify subplot
                subplot(2,2,i);

                %Calculate interpolated value 
                d_test = mean(params.break.funcs_d_range(i,:));
                eps_test = mean(params.break.funcs_eps_range(i,:));
                paramsin.turb.eps(:) = eps_test;
                [beta_eddy_test, beta_ratio_test, b_eddy_test] = BreakageInterpolate(d_test, eps_test, 1, params);

                %Calculate actual value for comparison
                if strcmp(params.break.model, 'Luo_Svendson_1996')
                    [b_eddy, beta, int_ratio] = BreakageLuoSvendson(1, 1, d_test, 1, lambda_min, N_lambdas, paramsin);
                elseif strcmp(params.break.model, 'Wang_2005')
                    [b_eddy, beta, int_ratio] = BreakageEddyAlt(1, 1, d_test, 1, lambda_min, N_lambdas, paramsin);
                else
                    error('Invalid breakage model. Options are "Wang_2005", and "Luo_Svendson_1996".');
                end

                %Compare the values
                b_eddy_err = abs(b_eddy - b_eddy_test)./b_eddy;
                beta_err = mean(abs(beta - beta_eddy_test)./beta);

                %Plot the betas
                plot(params.fvs_norm_all, beta, 'b-', 'LineWidth', 1.5); hold on;
                plot(params.fvs_norm_all, beta_eddy_test, 'r.', 'MarkerSize', 18); hold on;
                xlabel('Breakage Fraction'); ylabel('BSD'); 
                title(sprintf('d = %.4f m, eps = %.2f m^2/s^3', d_test, eps_test));
                grid on; grid minor; axis square;
                ylim([0, 20]);
                legend('Actual', 'Interpolated');
                
                x = 1;
            end


        % catch
        %     error('Issue with interpolation functions. Review code.');
        % end

end