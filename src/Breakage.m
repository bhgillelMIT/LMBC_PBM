function [b_src, b_snk, b_mats] = Breakage(y, params)

    %Make parameters global
    global params

    %Debug message
    if params.debug
        t_start = cputime;
        fprintf('--Breakage Start\n');
    end

    %Settings
    params.delta = 1E-3; %Value of 0.01 is used in Wnag et al (2003)


    %debug dettings
    if params.debug
        %params.turb.eps = 1 .* ones(size(params.turb.eps));
    end

    %Pull values and allocate output vectors
    b_src = zeros(size(y)); b_snk = b_src;
    b_mats = cell(params.Nz, 1);
    

    %Calculate kolmogorov length scale
    lambda_komogorov = ((params.nus.^3)./params.turb.eps).^0.25; %m
    lambda_min = 31.4 * lambda_komogorov; %Minimum eddy diameter to break a bubble - integrate up to bubble diameter 

    %Determien iz inds to use
    if params.sol.sep_layer & strcmp(params.sol.type, 'segregated')
        iz_inds = params.iz;
    else
        iz_inds = 1:params.Nz;
    end


    %Iterate through all spatial cells
    for iz = iz_inds

        %Pull numeric densities in this spatial celll
        z = params.zms(iz);
        cellinds = ((iz-1).*params.Nms + 1):1:(iz.*params.Nms);
        Ns = y();
        Ns_cell = Ns(cellinds);

        %Pull local fluid properties
        uL = 0;

        %Calculate critical breakup diameter
        d_crit = CalcCritDiameter(iz, params); %m - critical diameter for instability breakage

        %Create storage matrix
        bs = zeros(params.Nms, params.Nms); %Row = bubble breaking, column = bubble receiving
        b_mat = zeros(params.Nms);


        %VARY RESOLUTION OF fvs norm based on the number of size groups
        %below it to reduce computational cost

        %Iterate through size groups
        betas = zeros(params.Nms-1, length(params.fvs_norm_all));
        for im = 2:params.Nms %params.Nms %Neglect first size group, since no where for it to break

            %Pull index
            mind = cellinds(im);

            %Pull bubble size
            mi = params.mms(im); %Current represenative mass
            ni = params.nms(im);
            V = (ni .* (1 + params.X_mu(iz, im)) .* params.R .* params.T_mu(iz,im))./params.p_func(z);
            d = (6.*V/pi).^(1/3);
            %u = params.ubs.funcs{iz}(d); %params.uzs(cellinds(im));
            
            %Eddy Breakage
            if params.break.eddy
                switch params.break.model
                    case 'Luo_Svendson_1996'
                        [b_eddy, beta_eddy] = BreakageLuoSvendson(iz, im, d, Ns_cell, lambda_min, params.break.N_lambdas, params);
                        beta_eddy = beta_eddy./2;
                        %[beta_wang, beta_ratio, ~] = BreakageInterpolate(d, params.turb.eps(iz), iz, params);
                        %beta_eddy = beta_wang./trapz(params.fvs_norm_all, beta_wang);

                    case 'Wang_2005'


                        if params.break.interp %Solve using interpolation
                            if any(~isreal(params.turb.eps(iz)))
                                x = 1;
                            end
        
                             ts = cputime; 
        
        
                            [beta_eddy, beta_ratio, b_eddy] = BreakageInterpolate(d, params.turb.eps(iz), iz, params);
        
                            % %if params.sol.solve_details
                            % beta_ratio = params.break.beta_ratio{iz}(d, params.turb.eps(iz));
                            %     %b_eddy = params.break.b_eddy{iz}(d, params.turb.eps(iz));
                            % 
                            % tm  = cputime;
                            % tmd = tm - ts;
        
                            
                            % 
                            % beta_eddy = params.break.beta{iz}(d, params.turb.eps(iz), params.fvs_norm_all);
                            % beta_eddy = beta_eddy(:);
                            if any(isnan(beta_eddy))
                                beta_eddy = zeros(size(beta_eddy));
                                beta_ratio = 0;
                            end
        
                            %Log result 
                            params.break.beta_eddy_last = beta_eddy;
                            params.break.beta_ratio_last = beta_ratio;
                            % 
                            % 
                            %   te = cputime;
                            % td = te - ts;
                           % else
                            %    beta_eddy = params.break.beta_eddy_last;
                            %    beta_ratio = params.break.beta_ratio_last;
                           % end
        
                            % if td > 1E-8
                            %     x =1;
                            % end

                            b_eddy = BreakageEddySimple(iz, im, d, Ns_cell, beta_ratio, beta_eddy, params);


                        else %Solve in real time
                            [b_eddy, beta_eddy] = BreakageEddyAlt(iz, im, d, Ns_cell, lambda_min, params.N_lambdas, params); %b_eddy = breakage rate of the size; beta_eddy = distribution of bubble sizes
                        end
                    case 'Uniform_Binary'
                        [b_eddy, beta_eddy] = BreakageEddyUniformBinary(iz, im, d, Ns_cell, params);

                    otherwise
                        error('Invalid Breakage Model.')

                end
            end
            % 
            % %
            % % 
            % if strcmp(params.break.model, 'Uniform_Binary')
            %     beta_eddy_orig = beta_eddy;
            %     beta_eddy = ones(size(beta_eddy));
            % end

            %Inertia Breakage - Bubble too large
            if params.break.surf
                if d > d_crit
                    [b_surf, beta_surf] = BreakageSurf(d, d_crit, params);
                else
                    b_surf = 0;
                    beta_surf = 0;
                end
            else
                b_surf = 0;
                beta_surf = 0;
            end

            %Log BSD
            betas(im-1, :) = beta_eddy;

            if b_eddy > 0
                x = 1;
            end

            %SOURCE TERM --------------------------------------------------
            b_eddy = b_eddy * Ns_cell(im);
            b_surf = b_surf .* Ns_cell(im);
            b_total = b_eddy + b_surf;

            if b_surf > 0
                x = 1;
            end

            %Only consider if breakage rate is above a threshold
            if b_total > 1E-6 

                %Calculate distribution of products
                ms_norm = mi * params.fvs_norm_all;
                ms_max = max(ms_norm);
                    %ds_norm = (d.^3 .* params.fvs_norm).^(1/3); %Can calculate and list as a matrix with N_fvs columns, and Nms rows 
                zetas = zeros(1,im);

                %Create interpolation function
                %beta_func = @(m) (1/ms_max) .* interp1(ms_norm, beta_eddy, m); %Normalized to mass - CHECK THIS TO SOLVE MASS BALANCE ISSUE
                %beta_func = @(m) interp1(ms_norm, beta_eddy, m);
                beta_func = @(m) (1/mi) .* interp1(ms_norm, beta_eddy, m);
                beta_func_int = integral(beta_func, 0, mi);
                beta_func = @(m) 2/beta_func_int .* beta_func(m);
    
                %Iterate through smaller brackets
                for is = 1:im %Only consider sizes smaller than 


                    %SEARCH HERE FOR MISTAKE
    
                    %Pull cellind
                    sind = cellinds(is); %source index
    
                    


                    %Specify representative masses of adjacent cells (with
                    %   which the new population is shared)
                    if is > 1
                        m_low = params.mms(is-1);
                    else
                        m_low = 0;
                    end
                    if is < params.Nms
                        m_hig = params.mms(is+1);
                    end
                    m_mid = params.mms(is);
                    
                    %Calculate integral
                    %if im == 1 && is == 1
    
                    if is == 1 %Only preserve mass, not numbers
                        int_func_up = @(m) (m_hig - m)./(m_hig - m_mid) .* beta_func(m);
                        int_func_down = @(m) (m - 0)./(m_mid) .* beta_func(m); %Will conserve mass, but not number of bubbles since there is no lower bracket
                        zetas(is) = integral(int_func_up, m_mid, m_hig) + integral(int_func_down, 0, m_mid);
                    elseif is == params.Nms
                        int_func_down =  @(m) (m - m_low)./(m_mid - m_low) .* beta_func(m);
                        zetas(is) = integral(int_func_down, m_low, m_mid); %Don't consider up direction since bubbles cannot be larger;
                    elseif is == im 
                        int_func_down =  @(m) (m - m_low)./(m_mid - m_low) .* beta_func(m);
                        zetas(is) = integral(int_func_down, m_low, m_mid); %Don't consider up direction since bubbles cannot be larger;
                    else
                        int_func_up = @(m) (m_hig - m)./(m_hig - m_mid) .* beta_func(m);
                        int_func_down = @(m) (m - m_low)./(m_mid - m_low) .* beta_func(m);
                        zetas(is) = integral(int_func_up, params.mms(is), params.mms(is+1)) + integral(int_func_down, params.mms(is-1), params.mms(is));
                    end
    
                    %Add to source term
                    b_src_s = zetas(is) .* b_eddy + b_surf .* 1/(im-1);
                    b_src(sind) = b_src(sind) + b_src_s; % 1/(im-1) * %UPDATE UPDATE UPDATE
                    
                    %Store 
                    b_mat(im, is) = b_mat(im, is) + b_src_s;
                    

                end


            elseif strcmp(params.break.model, 'Equal_Binary') && b_total > 1E-6
                
                %Identify bins to split the new bubbles between
                m_break = mi/2;
                ind_lower = max(find(params.mms < m_break));
                ind_upper = min(find(params.mms >= m_break));

                %Handle case for smallest bin
                if isempty(ind_lower)
                    m_ratio = params.mms(1)/m_break;
                    eta_lower = 1/m_ratio;
                    eta_upper = 0;

                else
                    eta_lower = (params.mms(ind_upper) - m_break)/(params.mms(ind_upper) - params.mms(ind_lower));
                    eta_upper = 1 - eta_lower;
         
                end

                %Allocate bubble
                b_src(ind_lower) = eta_lower * 2 * b_eddy;
                b_src(ind_upper) = eta_upper * 2 * b_eddy;

                
                x = 1;

                



            end

            %SINK TERM ----------------------------------------------------
            b_snk(mind) = b_total;

            
           

            % if sum(m_src) > 0
            %     b_src = (sum(m_snk)./sum(m_src)) .* b_src;
            % end

        end

         %Check mass conservation and rescale 
        if params.src.debug
            params.break.m_src(params.src.its) = sum(b_src .* [repmat(params.mms, 1, params.Nz)]); %params.mms_rep); %[repmat(params.mms, 1, params.Nz)]'); %Derived from the distribution equations
            params.break.m_snk(params.src.its) = sum(b_snk .* [repmat(params.mms, 1, params.Nz)]); %params.mms_rep); %[repmat(params.mms, 1, params.Nz)]'); %Calculated
                
            %Normalize distributions
            b_src = (params.break.m_snk(params.src.its)./(params.break.m_src(params.src.its) + 1E-16)) .* b_src;
        
        end

        %Debug plots
        if params.break.debug & false

            %Create initial figure
            figure('units', 'normalized', 'OuterPosition',  [0, 0, 1, 1]);
            linespecs = {'-', '--', ':', '-.'};
            it = 2;
            for i = 2:5; %params.Nms

                %Plot BSD
                %subplot(3,5,it);
                plot(params.fvs_norm_all, betas(i-1,:), 'k-', 'LineWidth', 1.5, 'LineStyle', linespecs{i-1}); hold on;
                grid on; grid minor; axis square;
                xlabel('fv - breakup vol. fraction'); ylabel('Normalized BSD');
                title(sprintf('BSD (im = %d)', i));
                it = it + 1;

                %Create new figure and reset iteration counter
                if it > 15
                    figure('units', 'normalized', 'OuterPosition',  [0, 0, 1, 1]);
                    it = 1;
                end
            end

            %Plots
            title('Normal Bubble Size Distribution');
            legend('d = 1.5 mm', 'd = 2 mm', 'd = 3 mm', 'd = 6 mm');
            set(gca, 'FontSize', 18);


        end

        %Store source matrix
        b_mats{iz} = b_mat;

    end

    %Debug message
    if params.debug
        t_end = cputime;
        t_req = t_end - t_start;
        fprintf('--Breakage End (t = %0.4f s)\n', t_req);
    end

    % if params.sol.solve_details
    %     params.sol.solve_details = false;
    % end


end



function [b_surf, fv] = BreakageSurf(d, d_crit, params)

    %Additional Parameters - pronably move outside loop
    

    %
    b_surf = params.break.b_star .* ((d - d_crit).^params.break.m_star)./((d - d_crit).^params.break.m_star + d_crit.^params.break.m_star);
    fv = 0.5; %0.5 = equal breakup

end


function [b_eddy, beta] = BreakageEddy(iz, im, d, Ns_cell, lambda_min, N_lambdas, params)

    %Debug cases
    if params.break.debug
        %d = 0.003;
        t_start.total = cputime;
        fprintf('-Eddy Breakage: iz = %i; im = %i; d = %0.5f m; m_i = %0.3e;\n', iz, im, d, params.mms(im));
        x = 1;
    end

    %Specify standard fv distribution to interpolate onto
    % fvs_norm = [0:params.dfv_cap:0.05, (0.05+params.dfv_surf):params.dfv_surf:0.5];

    %Create storage vectors
    Pbs_norm = params.Pbs_norm; %zeros(N_lambdas, length(params.fvs_norm));

    %Iterate through lambdas
    it_total = 0;
    if lambda_min(iz) < d %otherwise the bubble cannot be broken by eddies

        %
        lambdas = linspace(lambda_min(iz), d, N_lambdas);

        %Calculate cutoffs for each lambda
        u_bar_lambdas = sqrt(2) .* (params.turb.eps(iz) .* lambdas).^(1/3);
        e_bar_lambdas = (pi/6) .* lambdas.^3 .* params.rhos(iz) .* (u_bar_lambdas.^2)./2; %Might be able to move this calculation outside of ODE func - use interp
        e_cutoffs = 10 .* e_bar_lambdas;

        for il = 1:N_lambdas

            %Debug timing
            if params.debug
                t_start.setup = cputime;
            end
    
            %Pull current eddy diamter
            lambda = lambdas(il);
        
            lambda_rat =  lambda/d;
            %lambda = lambda_rat * d;
    
            %Evaluate transition break fraction 
            fvc = params.fvc_func(lambda_rat);
    
            %Specify values to integrate over
            fvs_surf = 0.5:(-params.dfv_surf):fvc; fvs_surf(end+1) = fvc;
            fvsm_surf = (fvs_surf(1:end-1) + fvs_surf(2:end))./2;
            if fvc/params.dfv_cap > params.N_fv_min
                fvs_cap = 1E-8:params.dfv_cap:fvc;
            else
                fvs_cap = linspace(1E-8,fvc,params.N_fv_min);
            end
            fvsm_cap = (fvs_cap(1:end-1) + fvs_cap(2:end))./2;
    
            %Calculate average velocity and energy of eddy
            %u_bar_lambda = sqrt(2) .* (params.turb.eps(iz) .* lambda).^(1/3);
            %e_bar_lambda = (pi/6) .* lambda.^3 .* params.rhos(iz) .* (u_bar_lambda.^2)./2;
            %e_cutoff = 10 * e_bar_lambda;

            %Calculate min and max breakup fraction
            cf_max = min([(2^(1/3) - 1), e_cutoffs(il)./(pi .* d.^2 .* params.sigmas(iz))]);
            fv_func = @(fv) fv.^(2/3) + (1-fv).^(2/3) - 1 - cf_max;
            fv_max = fzero(fv_func, [0, 0.5]);
            fv_max = round(fv_max,3);

            if fv_max > 0.49
                x = 1;
            end

            % fv_min = ((pi .* lambda.^3.*params.sigmas(iz))./(6 * e_cutoff * d)).^3;
            % cf_max = min([(2.^(1/3) - 1), e_cutoff./(pi .* d.^2 .* params.sigmas(iz))]);
            % fv_func = @(fv) fv.^(2/3) + (1-fv).^(2/3) - 1 - cf_max;
            % fv_max = fzero(fv_func, 0.5);
            % cf_half = 2^(1/3) - 1;
            % e_half = cf_half * pi * d^2 * params.sigmas(iz);
            
            %Calculate energy required to break the bubble into
            %each possible fraction
            cfs = (fvs_surf.^(2/3) + (1-fvs_surf).^(2/3) - 1);
            e_is = cfs .* pi .* d.^2 .* params.sigmas(iz); %Energy required to break - should be calculated outside of loop
            e_ihs = (e_is(1:end-1) + e_is(2:end))./2;

            %Calculate initial value at 0.5
            if e_is(1) < e_cutoffs(il)
                fv = 0.5;
                e_min = (fv.^(2/3) + (1-fv).^(2/3) - 1) .* pi .* d.^2 .* params.sigmas(iz);
                e_max = e_cutoffs(il);
                fv_min = ((pi .* lambda.^3 .* params.sigmas(iz))./(6 .* e_ihs(1) .* d)).^3;
                Pbe = 1./(fv - fv_min);
                Pe = (1./e_bar_lambdas(il)) .* exp(-e_ihs(1)/e_bar_lambdas(il)); %
                Pb_surf_init = Pbe .* Pe .* (e_is(1) - e_is(2));


                P_func = @(e) 1./(fv - fv_min) .* (1./e_bar_lambdas(il)) .* exp(-e/e_bar_lambdas(il));

                Pb_surf_init = integral(P_func, e_is(1), e_max);

            else
                Pb_surf_init = 0;
            end


             %To calculate initial Pb at fv = 0.5 for the surface controlled region, need to consider cases fv_max = 0.5. If fv_max < 0.5, then it starts at zero as it currently does, otherwise the baseline is higher.

            %Debug timing
            if params.debug
                t_end = cputime;
                t_req.setup = t_end - t_start.setup;
                fprintf('--Eddy Breakage | Setup time   = %0.6f s\n', t_req.setup);
                t_start.surf = cputime;
            end

            %Iterate through possible breakup fractions
            Pb_surf = zeros(1, length(fvsm_surf));
            can_break_surf = e_ihs <  e_cutoffs(il); %Logical stating if any eddies have the energy to break the bubble into the defined fractions
    
            %Surface energy region ----------------------------------------
            if any(can_break_surf) %Check if no eddies can break this bubble
                
                

                %Iterate only through cases which can cause
                %breakage
                break_inds = find(can_break_surf);
                for iv = break_inds

                    %Calculate probability of eddy having the required kinetic energy
                    Pe = (1./e_bar_lambdas(il)) .* exp(-e_ihs(iv)/e_bar_lambdas(il)); %probability of an eddy with this amount of energy

                        %Vectorize and do outside loop
    
                    fv = fvs_surf(iv);
                    %fvh = (fvs_surf(iv) + fvs_surf(iv+1))./2;
                    %cf = (fv.^(2/3) + (1-fv).^(2/3) - 1);
                    %efv = cf * pi * d^2 * params.sigmas(iz);
                    
    
                    %Calculate minimum energy and minimum breakup fraction
                    e_min = (fv.^(2/3) + (1-fv).^(2/3) - 1) .* pi .* d.^2 .* params.sigmas(iz);
                    fv_ih_min = ((pi .* lambda.^3 .* params.sigmas(iz))./(6 .* e_ihs(iv) .* d)).^3;
    
                    %Calculate breakage probability
                    if (fvsm_surf(iv) - fv_ih_min) >= params.break.delta
                        Pbe = 1/(fvs_surf(iv) - fv_ih_min);
                    else
                        Pbe = 0;
                    end
    
                    %Calculate total probability
                    if iv == 1
                        Pb_surf(iv) = (e_is(iv) - e_is(iv+1)) .* Pbe .* Pe + Pb_surf_init;
                    else
                        %Pb_surf(iv) = (e_is(iv) - e_is(iv+1)) .* Pbe .* Pe + Pb_surf(iv-1);

                        %Alternative calculation
                        P_func = @(e) 1./(fv - fv_ih_min) .* (1./e_bar_lambdas(il)) .* exp(-e/e_bar_lambdas(il));
                        Pb_surf(iv) = integral(P_func, e_ihs(iv), e_cutoffs(il)) + Pb_surf(iv-1);
                    end  

                    

                    %Log number of iterations
                    it_total = it_total + 1;
    
                end
    
                %Fill rest of 
            end

            %Debug timing
            if params.debug
                t_end = cputime;
                t_req.surf= t_end - t_start.surf;
                fprintf('--Eddy Breakage | Surf. time   = %0.6f s\n', t_req.surf);
                t_start.cap = cputime;
            end
            
    
    
    
            %Pressure controlled region --------------------------
            e_is = (pi .* lambda.^3 .* params.sigmas(iz))./(6 .* fvs_cap.^(1/3) .* d);
            e_ihs = (e_is(1:end-1) + e_is(2:end))./2;
    
            %Iterate through possible breakup fractions
            Pb_cap = zeros(1, length(fvsm_cap));
            can_break_cap = e_ihs < e_cutoffs(il); 
            if any(can_break_cap)
                break_inds = find(can_break_cap);
                for iv = break_inds
    
                    %Calculate Pbe - Pb for the specific energy
                    fv_ih_max_func = @(fv) e_ihs(iv) - (fv.^(2/3) + (1 - fv).^(2/3) - 1) .* pi .* d.^2 .* params.sigmas(iz);
                    try
                        fv_ih_max = fzero(fv_ih_max_func, [0, 0.5]);
                        %fv_ih_max = params.fvih_max_funcs{iz}(d, e_ihs(iv)); %fzero(fv_ih_max_func, 0.05);
                        if fv_ih_max > 0 && (fv_ih_max - fvsm_cap(iv)) >= params.delta
                            Pbe = 1/(fv_ih_max - fvsm_cap(iv));
                        else
                            Pbe = 0;
                        end
                    catch
                        Pbe = 0;
                    end
    
                    %Calculate 
                    Pe = (1./e_bar_lambdas(il)) .* exp(-e_ihs(iv)/e_bar_lambdas(il));
    
                     %Calculate total probability
                    if iv == 1
                        Pb_cap(iv) = (e_is(iv) - e_is(iv+1)) .* Pbe .* Pe;
                    else
                        Pb_cap(iv) = (e_is(iv) - e_is(iv+1)) .* Pbe .* Pe + Pb_cap(iv-1);
                    end  
    

                    if any(Pb_cap < 0)
                        x = 1;
                    end

                    %Log number of iterations
                    it_total = it_total + 1;
    
                end
    
            end

            %Debug timing
            if params.debug
                t_end = cputime;
                t_req.cap= t_end - t_start.cap;
                fprintf('--Eddy Breakage | Cap. time    = %0.6f s\n', t_req.cap);
                t_start.interp = cputime;
            end
            
            
    
            %Interpolate results
            fvs = [fvsm_cap, fliplr(fvsm_surf)];
            [fvs, unique_inds] = unique(fvs);
            Pbs = [Pb_cap, fliplr(Pb_surf)];
            Pbs = Pbs(unique_inds);
            Pbs_norm(il,:) = interp1(fvs, Pbs, params.fvs_norm, 'linear', 'extrap');
            Pbs_norm(il,1) = 0;
            % 
            %Debug analysis
            if params.debug

                %Plot
                figure();
                plot(params.fvs_norm, Pbs_norm(il,:), 'k-', 'LineWidth', 2);
                grid on; grid minor; axis square;

                if lambda_rat > 0.69
                    x =1;
                end
                close 

            end

            %Debug timing
            if params.debug
                t_end = cputime;
                t_req.interp = t_end - t_start.interp;
                fprintf('--Eddy Breakage | Interp. time = %0.6f s\n', t_req.interp);
            end

            
    
            %Calculate complete integral
            %int_lambda(il) = sum(Pbs_norm(il,:)) .* ((lambdas(il) + d).^2)./(lambdas(il).^(11/3)); 
    

    
            %  %Debug updates
            % if params.debug
            %     fprintf('im = %d; d = %0.4f mm; iv = %d; Pb = %0.4e;\n', im, d, iv, Pb)
            % end
        end

    else

        %Probability is zero - lambda is too large, will just transport
        %bubble 
        error('Finish this.');
        
    end


    %Evaluate outer integral
    b_fvd = params.bfd_zero;
    for iv = 2:length(params.fvs_norm)
        f_lambda_fv = Pbs_norm(:,iv)' .* ((lambdas + d).^2)./(lambdas.^(11/3));
        int_lambda_fv = trapz(lambdas, f_lambda_fv);
        b_fvd(iv) = 0.923 .* (1 - params.alpha_g(iz)) .* Ns_cell(im) .* params.turb.eps(iz).^(1/3) .* int_lambda_fv;

        %Log number of iterations
        it_total = it_total + 1;
    end

    %Calculate overall rate
    b_eddy = trapz(params.fvs_norm, b_fvd);

    %Determine normalized daughter bubble size distribution
    beta = (b_fvd)./(b_eddy);
    beta = [beta, fliplr(beta(1:end-1))];

    %
    if any(isnan(beta))
        beta = zeros(size(beta));
    end


    %Debug plot
    if params.break.debug
        t_end = cputime;
        t_req = t_end - t_start.total;
        figure();
        plot(params.fvs_norm_all, beta);
        close
        fprintf('-- Eddy Breakage time = %0.4f; Total its = %d; \n', t_req, it_total)
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

    if beta_ratio > 0
        x = 1;
    end

    %Normalize beta_eddy
    beta_eddy_int = trapz(params.fvs_norm_all, beta_eddy);
    beta_eddy = beta_eddy./beta_eddy_int;

    %Track iterations
    it_total = 1;

    %Calculate integral
    b_fvd = params.break.bfd_zero;
    int_lambda_fvs = beta_ratio .* beta_eddy(1:length(params.fvs_norm));
    b_fvd = 0.923 .* (1 - params.alpha_g(iz))  .* params.turb.eps(iz).^(1/3)  .* int_lambda_fvs; %.* Ns_cell(im)

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

    if any(b_fvd > 0)
        x = 1;
    end

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

    if b_eddy > 0
        x = 1;
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



function [b_eddy, beta] = BreakageEddyUniformBinary(iz, im, d, Ns_cell, params)

    %Pull constant breakage rate
    b0 = params.break.constant_rate;

    %Calculate breakage rate
    b_eddy = b0 * (params.mms(im).^2) .* Ns_cell(im);

    %Define beta
    %beta = 2./params.mms(im) .* ones(size(params.fvs_norm_all));
    beta = ones(size(params.fvs_norm_all));
    % beta = zeros(size(params.fvs_norm_all));
    % ind_mid = (length(params.fvs_norm_all) + 1)/2;
    % beta(ind_mid) = 1./(params.fvs_norm_all(ind_mid) - params.fvs_norm_all(ind_mid-1));



end