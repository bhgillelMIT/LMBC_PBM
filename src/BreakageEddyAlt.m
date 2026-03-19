function [b_eddy, beta, int_ratio] = BreakageEddyAlt(iz, im, d, Ns_cell, lambda_min, N_lambdas, params)

    %Debug cases
    if params.debug
        t_start.total = cputime;
        fprintf('-Eddy Breakage (Wang): iz = %i; im = %i; d = %0.5f m; m_i = %0.3e;\n', iz, im, d, params.mms(im));
        x = 1;
    end

    %Create storage vectors
    Pbs_norm = params.Pbs_norm; %zeros(N_lambdas, length(params.fvs_norm));
    its_total = 0;
    if lambda_min(iz) < d %otherwise the bubble cannot be broken by eddies

        %List lambdas 
        lambdas = linspace(lambda_min(iz), d, N_lambdas);


        e_cap_min = (pi .* lambdas.^3 .* params.sigmas(iz))./(6 .* 0.5.^(1/3) .* d);
        e_surf_min = (2^(1/3) - 1) .* pi .* d.^2 .* params.sigmas(iz);
        

        %Determine mean energy of eddies with each lambda
        u_bar_lambdas = sqrt(2) .* (params.turb.eps(iz) .* lambdas).^(1/3);
        e_bar_lambdas = (pi/6) .* lambdas.^3 .* params.rhos(iz) .* (u_bar_lambdas.^2)./2; %Might be able to move this calculation outside of ODE func - use interp
        e_cutoffs = 10 .* e_bar_lambdas;

        %Iterate through fvs
        fvhs = (params.fvs_norm(1:end-1) + params.fvs_norm(2:end))./2;
        for iv = 1:length(fvhs)
            fv = fvhs(iv);
            
            %Determine minimum lambda
            e_cap_min = (pi .* lambdas.^3 .* params.sigmas(iz))./(6 .* fv.^(1/3) .* d);
            cf = (fv.^(2/3) + (1-fv).^(2/3) - 1);
            e_surf_min = cf .* pi .* d.^2 .* params.sigmas(iz);
            
            Pbs = zeros(size(lambdas));
            for il = 1:length(lambdas)
                
                %Pull lambda
                lambda = lambdas(il);

                %Calculate cticial volume fraction
                lambda_rat = lambda/d;
                fvc = params.fvc_func(lambda_rat);

                %Calculate energy required
                cf = (fv.^(2/3) + (1-fv).^(2/3) - 1);
                e_surf = cf .* pi .* d.^2 .* params.sigmas(iz);
                e_cap  = (pi .* lambda.^3 .* params.sigmas(iz))./(6 .* fv.^(1/3) .* d);
                e_rat = e_surf/e_cap;

                %Calculate integral
                P_func = @(e) Pb_func(e, lambda, params, fv, d, iz, e_bar_lambdas(il));
                Pbs_norm(il, iv) = integral(P_func, 1E-12, e_cutoffs(il));

                %Log iterations
                its_total = its_total + 1;
    
    
            end



        end

        %Evaluate outer integral
        if any(Pbs_norm(:) > 0)
            b_fvd = params.break.bfd_zero;
            ints_fvd = params.break.bfd_zero;
            for iv = 2:length(params.fvs_norm)
                f_lambda_fv = Pbs_norm(:,iv)' .* ((lambdas + d).^2)./(lambdas.^(11/3));
                int_lambda_fv = trapz(lambdas, f_lambda_fv);
                ints_fvd(iv) = int_lambda_fv;
                b_fvd(iv) = 0.923 .* (1 - params.alpha_g(iz)) .* Ns_cell(im) .* params.turb.eps(iz).^(1/3) .* int_lambda_fv;
                
                %Log number of iterations
                its_total = its_total + 1;
            end

        else
            b_fvd = params.break.bfd_zero;
            ints_fvd = b_fvd;
        end
        
        %Calculate overall rate
        b_eddy = trapz(params.fvs_norm, b_fvd);

        %Determine normalized daughter bubble size distribution
        beta = (b_fvd)./(b_eddy);

        %Handle NaNs
        if any(isnan(beta))
            beta = zeros(size(beta));
        end

        %Extrapolate to fill in center value
        try
            interp_func = @(fv) interp1(params.fvs_norm(1:end-1), beta(1:end-1), fv, 'spline');
            beta(end) = interp_func(0.5);
            beta = [beta, fliplr(beta(1:end-1))];
        catch
            x = 1;
        end

        %Handle NaNs
        if any(isnan(beta))
            beta = zeros(size(beta));
        end

        %Calculate ratio of ints_fvd to beta
        int_ratio = ints_fvd./beta(1:length(ints_fvd));
        int_ratio = mode(int_ratio);
        if isnan(int_ratio)
            int_ratio = 0;
        end

        if b_eddy > 0
            x = 1;
        end
    
        %Test that the int_ratio converts back to the breakage rate
        int_lambda_fvs = int_ratio .* beta(1:length(params.fvs_norm));
    
        b_fvd_test = 0.923 .* (1 - params.alpha_g(iz)) .* Ns_cell(im) .* params.turb.eps(iz).^(1/3)  .* int_lambda_fvs;
        b_eddy_test = trapz(params.fvs_norm, b_fvd_test);

        %Check b_eddy_test  
        b_eddy_test_diff = abs(b_eddy_test - b_eddy);
        if b_eddy_test_diff > (0.0001 * b_eddy)
            x = 1;
        end

    else
        b_eddy = 0;
        int_ratio = 0;
        beta = ones(size(params.fvs_norm_all));

    end



    %Gap at center
    Pbs_norm(:,end) = Pbs_norm(:,end-1);


    

    
    

   

    

    

    %Debug plot
    if params.debug
        t_end = cputime;
        t_req = t_end - t_start.total;
        figure();
        plot(params.fvs_norm_all, beta);
        close
        fprintf('-- Eddy Breakage time = %0.4f; Total its = %d; \n', t_req, its_total)
    end


   
end




function Pb = Pb_func(e, lambda, params, fv, d, iz, e_bar_lambda)

    %Calculate minimum size
    fv_min = ((pi .* lambda.^3 .* params.sigmas(iz))./(6 .* e.* d)).^3;
    cf_max_max = (2^(1/3) - 1) .* ones(size(e));
    cf_max = min([cf_max_max; e./(pi .* d.^2 .* params.sigmas(iz))]);

    %Feasible inds
    test_inds = find(fv_min < 1);
    
    %Iterate through
    Pb = zeros(size(e));

    

    for ie = test_inds %1:length(e)

        %Calculate fv_max
        if cf_max(ie) == cf_max_max(ie)
            fv_max = 0.5;
        else
            fv_max_func = @(fv) cf_max(ie) - (fv.^(2/3) + (1 - fv).^(2/3) - 1);
            try
                fv_max = fzero(fv_max_func, [0,0.5]); %params.fvih_max_funcs{iz}(d, e(ie)); %fzero(fv_max_func, [0,0.5]);
            catch
                fv_max = 0;
            end
        end


        %CONSIDER IF 

        if (fv_max - fv_min(ie)) > params.break.delta && (fv < fv_max) && (fv > fv_min(ie))
    
            Pb(ie) = 1./(fv_max - fv_min(ie)) .* (1./e_bar_lambda) .* exp(-e(ie)/e_bar_lambda);
    
        else
            Pb(ie) = 0;
        end
    end

    if any(Pb > 0)
        x =1;
    end

end





    % %Iterate through lambdas
    % it_total = 0;
    % if lambda_min(iz) < d %otherwise the bubble cannot be broken by eddies
    % 
    %     %Create discrete lambdas for integration
    %     lambdas = linspace(lambda_min(iz), d, N_lambdas);
    % 
    %     %Calculate cutoffs for each lambda
    %     u_bar_lambdas = sqrt(2) .* (params.turb.eps(iz) .* lambdas).^(1/3);
    %     e_bar_lambdas = (pi/6) .* lambdas.^3 .* params.rhos(iz) .* (u_bar_lambdas.^2)./2; %Might be able to move this calculation outside of ODE func - use interp
    %     e_cutoffs = 10 .* e_bar_lambdas;
    % 
    %     for il = 1:N_lambdas
    % 
    %         %Debug timing
    %         if params.debug
    %             t_start.setup = cputime;
    %         end
    % 
    %         %Pull current eddy diamter
    %         lambda = lambdas(il);
    % 
    %         lambda_rat =  lambda/d;
    %         %lambda = lambda_rat * d;
    % 
    %         %Evaluate transition break fraction 
    %         fvc = params.fvc_func(lambda_rat);
    % 
    %         %Specify values to integrate over
    %         fvs_surf = 0.5:(-params.dfv_surf):fvc; fvs_surf(end+1) = fvc;
    %         fvsm_surf = (fvs_surf(1:end-1) + fvs_surf(2:end))./2;
    %         if fvc/params.dfv_cap > params.N_fv_min
    %             fvs_cap = 1E-8:params.dfv_cap:fvc;
    %         else
    %             fvs_cap = linspace(1E-8,fvc,params.N_fv_min);
    %         end
    %         fvsm_cap = (fvs_cap(1:end-1) + fvs_cap(2:end))./2;
    % 
    %         %Calculate average velocity and energy of eddy
    %         %u_bar_lambda = sqrt(2) .* (params.turb.eps(iz) .* lambda).^(1/3);
    %         %e_bar_lambda = (pi/6) .* lambda.^3 .* params.rhos(iz) .* (u_bar_lambda.^2)./2;
    %         %e_cutoff = 10 * e_bar_lambda;
    % 
    %         %Calculate min and max breakup fraction
    %         cf_max = min([(2^(1/3) - 1), e_cutoffs(il)./(pi .* d.^2 .* params.sigmas(iz))]);
    %         fv_func = @(fv) fv.^(2/3) + (1-fv).^(2/3) - 1 - cf_max;
    %         fv_max = fzero(fv_func, [0, 0.5]);
    %         fv_max = round(fv_max,3);
    % 
    %         if fv_max > 0.49
    %             x = 1;
    %         end
    % 
    %         % fv_min = ((pi .* lambda.^3.*params.sigmas(iz))./(6 * e_cutoff * d)).^3;
    %         % cf_max = min([(2.^(1/3) - 1), e_cutoff./(pi .* d.^2 .* params.sigmas(iz))]);
    %         % fv_func = @(fv) fv.^(2/3) + (1-fv).^(2/3) - 1 - cf_max;
    %         % fv_max = fzero(fv_func, 0.5);
    %         % cf_half = 2^(1/3) - 1;
    %         % e_half = cf_half * pi * d^2 * params.sigmas(iz);
    % 
    %         %Calculate energy required to break the bubble into
    %         %each possible fraction
    %         cfs = (fvs_surf.^(2/3) + (1-fvs_surf).^(2/3) - 1);
    %         e_is = cfs .* pi .* d.^2 .* params.sigmas(iz); %Energy required to break - should be calculated outside of loop
    %         e_ihs = (e_is(1:end-1) + e_is(2:end))./2;
    % 
    %         %Calculate initial value at 0.5
    %         if e_is(1) < e_cutoffs(il)
    %             fv = 0.5;
    %             e_min = (fv.^(2/3) + (1-fv).^(2/3) - 1) .* pi .* d.^2 .* params.sigmas(iz);
    %             e_max = e_cutoffs(il);
    %             fv_min = ((pi .* lambda.^3 .* params.sigmas(iz))./(6 .* e_ihs(1) .* d)).^3;
    %             Pbe = 1./(fv - fv_min);
    %             Pe = (1./e_bar_lambdas(il)) .* exp(-e_ihs(1)/e_bar_lambdas(il)); %
    %             Pb_surf_init = Pbe .* Pe .* (e_is(1) - e_is(2));
    % 
    % 
    %             P_func = @(e) 1./(fv - fv_min) .* (1./e_bar_lambdas(il)) .* exp(-e/e_bar_lambdas(il));
    % 
    %             Pb_surf_init = integral(P_func, e_is(1), e_max);
    % 
    %         else
    %             Pb_surf_init = 0;
    %         end
    % 
    % 
    %          %To calculate initial Pb at fv = 0.5 for the surface controlled region, need to consider cases fv_max = 0.5. If fv_max < 0.5, then it starts at zero as it currently does, otherwise the baseline is higher.
    % 
    %         %Debug timing
    %         if params.debug
    %             t_end = cputime;
    %             t_req.setup = t_end - t_start.setup;
    %             fprintf('--Eddy Breakage | Setup time   = %0.6f s\n', t_req.setup);
    %             t_start.surf = cputime;
    %         end
    % 
    %         %Iterate through possible breakup fractions
    %         Pb_surf = zeros(1, length(fvsm_surf));
    %         can_break_surf = e_ihs <  e_cutoffs(il); %Logical stating if any eddies have the energy to break the bubble into the defined fractions
    % 
    %         %Surface energy region ----------------------------------------
    %         if any(can_break_surf) %Check if no eddies can break this bubble
    % 
    % 
    % 
    %             %Iterate only through cases which can cause
    %             %breakage
    %             break_inds = find(can_break_surf);
    %             for iv = break_inds
    % 
    %                 %Calculate probability of eddy having the required kinetic energy
    %                 Pe = (1./e_bar_lambdas(il)) .* exp(-e_ihs(iv)/e_bar_lambdas(il)); %probability of an eddy with this amount of energy
    % 
    %                     %Vectorize and do outside loop
    % 
    %                 fv = fvs_surf(iv);
    %                 %fvh = (fvs_surf(iv) + fvs_surf(iv+1))./2;
    %                 %cf = (fv.^(2/3) + (1-fv).^(2/3) - 1);
    %                 %efv = cf * pi * d^2 * params.sigmas(iz);
    % 
    % 
    %                 %Calculate minimum energy and minimum breakup fraction
    %                 e_min = (fv.^(2/3) + (1-fv).^(2/3) - 1) .* pi .* d.^2 .* params.sigmas(iz);
    %                 fv_ih_min = ((pi .* lambda.^3 .* params.sigmas(iz))./(6 .* e_ihs(iv) .* d)).^3;
    % 
    %                 %Calculate breakage probability
    %                 if (fvsm_surf(iv) - fv_ih_min) >= params.delta_breakage
    %                     Pbe = 1/(fvs_surf(iv) - fv_ih_min);
    %                 else
    %                     Pbe = 0;
    %                 end
    % 
    %                 %Calculate total probability
    %                 if iv == 1
    %                     Pb_surf(iv) = (e_is(iv) - e_is(iv+1)) .* Pbe .* Pe + Pb_surf_init;
    %                 else
    %                     %Pb_surf(iv) = (e_is(iv) - e_is(iv+1)) .* Pbe .* Pe + Pb_surf(iv-1);
    % 
    %                     %Alternative calculation
    %                     P_func = @(e) 1./(fv - fv_ih_min) .* (1./e_bar_lambdas(il)) .* exp(-e/e_bar_lambdas(il));
    %                     Pb_surf(iv) = integral(P_func, e_ihs(iv), e_cutoffs(il)) + Pb_surf(iv-1);
    %                 end  
    % 
    % 
    % 
    %                 %Log number of iterations
    %                 it_total = it_total + 1;
    % 
    %             end
    % 
    %             %Fill rest of 
    %         end
    % 
    %         %Debug timing
    %         if params.debug
    %             t_end = cputime;
    %             t_req.surf= t_end - t_start.surf;
    %             fprintf('--Eddy Breakage | Surf. time   = %0.6f s\n', t_req.surf);
    %             t_start.cap = cputime;
    %         end
    % 
    % 
    % 
    % 
    %         %Pressure controlled region --------------------------
    %         e_is = (pi .* lambda.^3 .* params.sigmas(iz))./(6 .* fvs_cap.^(1/3) .* d);
    %         e_ihs = (e_is(1:end-1) + e_is(2:end))./2;
    % 
    %         %Iterate through possible breakup fractions
    %         Pb_cap = zeros(1, length(fvsm_cap));
    %         can_break_cap = e_ihs < e_cutoffs(il); 
    %         if any(can_break_cap)
    %             break_inds = find(can_break_cap);
    %             for iv = break_inds
    % 
    %                 %Calculate Pbe - Pb for the specific energy
    %                 fv_ih_max_func = @(fv) e_ihs(iv) - (fv.^(2/3) + (1 - fv).^(2/3) - 1) .* pi .* d.^2 .* params.sigmas(iz);
    %                 try
    %                     fv_ih_max = fzero(fv_ih_max_func, [0, 0.5]);
    %                     %fv_ih_max = params.fvih_max_funcs{iz}(d, e_ihs(iv)); %fzero(fv_ih_max_func, 0.05);
    %                     if fv_ih_max > 0 && (fv_ih_max - fvsm_cap(iv)) >= params.delta
    %                         Pbe = 1/(fv_ih_max - fvsm_cap(iv));
    %                     else
    %                         Pbe = 0;
    %                     end
    %                 catch
    %                     Pbe = 0;
    %                 end
    % 
    %                 %Calculate 
    %                 Pe = (1./e_bar_lambdas(il)) .* exp(-e_ihs(iv)/e_bar_lambdas(il));
    % 
    %                  %Calculate total probability
    %                 if iv == 1
    %                     Pb_cap(iv) = (e_is(iv) - e_is(iv+1)) .* Pbe .* Pe;
    %                 else
    %                     Pb_cap(iv) = (e_is(iv) - e_is(iv+1)) .* Pbe .* Pe + Pb_cap(iv-1);
    %                 end  
    % 
    % 
    %                 if any(Pb_cap < 0)
    %                     x = 1;
    %                 end
    % 
    %                 %Log number of iterations
    %                 it_total = it_total + 1;
    % 
    %             end
    % 
    %         end
    % 
    %         %Debug timing
    %         if params.debug
    %             t_end = cputime;
    %             t_req.cap= t_end - t_start.cap;
    %             fprintf('--Eddy Breakage | Cap. time    = %0.6f s\n', t_req.cap);
    %             t_start.interp = cputime;
    %         end
    % 
    % 
    % 
    %         %Interpolate results
    %         fvs = [fvsm_cap, fliplr(fvsm_surf)];
    %         [fvs, unique_inds] = unique(fvs);
    %         Pbs = [Pb_cap, fliplr(Pb_surf)];
    %         Pbs = Pbs(unique_inds);
    %         Pbs_norm(il,:) = interp1(fvs, Pbs, params.fvs_norm, 'linear', 'extrap');
    %         Pbs_norm(il,1) = 0;
    %         % 
    %         %Debug analysis
    %         if params.debug
    % 
    %             %Plot
    %             figure();
    %             plot(params.fvs_norm, Pbs_norm(il,:), 'k-', 'LineWidth', 2);
    %             grid on; grid minor; axis square;
    % 
    %             if lambda_rat > 0.69
    %                 x =1;
    %             end
    %             close 
    % 
    %         end
    % 
    %         %Debug timing
    %         if params.debug
    %             t_end = cputime;
    %             t_req.interp = t_end - t_start.interp;
    %             fprintf('--Eddy Breakage | Interp. time = %0.6f s\n', t_req.interp);
    %         end
    % 
    % 
    % 
    %         %Calculate complete integral
    %         %int_lambda(il) = sum(Pbs_norm(il,:)) .* ((lambdas(il) + d).^2)./(lambdas(il).^(11/3)); 
    % 
    % 
    % 
    %         %  %Debug updates
    %         % if params.debug
    %         %     fprintf('im = %d; d = %0.4f mm; iv = %d; Pb = %0.4e;\n', im, d, iv, Pb)
    %         % end
    %     end
    % 
    % else
    % 
    %     %Probability is zero - lambda is too large, will just transport
    %     %bubble 
    %     error('Finish this.');
    % 
    % end



                % 
                % 
                % %Calculate integrals in the surface region
                % if fvc < fv
                % 
                %     e_min = (fv.^(2/3) + (1-fv).^(2/3) - 1) .* pi .* d.^2 .* params.sigmas(iz);
                %     e_max = e_cutoffs(il);
                %     fv_min = ((pi .* lambda.^3 .* params.sigmas(iz))./(6 .* e_surf.* d)).^3;
                %     % Pbe = 1./(fv - fv_min);
                %     % Pe = (1./e_bar_lambdas(il)) .* exp(-e_ihs(1)/e_bar_lambdas(il)); %
                %     % Pb_surf_init = Pbe .* Pe .* (e_is(1) - e_is(2));
                %     P_func = @(e) 1./(fv - fv_min) .* (1./e_bar_lambdas(il)) .* exp(-e/e_bar_lambdas(il));
                % 
                %     Pb= integral(P_func, e_surf, e_max);
                % 
                % %Calculate ingrals in the capillary region
                % else
                %     cf_max = min([(2^(1/3) - 1), e_cap/(pi .* d.^2 .* params.sigmas(iz))]);
                %     fv_max = 
                % 
                % 
                %     try
                %         fv_ih_max = fzero(fv_ih_max_func, [0, 0.5]);
                % 
                % 
                % end
                % 