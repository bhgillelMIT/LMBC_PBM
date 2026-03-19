%This function implements the Luo and Svendson (1996) breakage model, which is based on the idea of a critical 
% eddy size that can cause breakage. The model calculates the breakage rate based on the turbulent energy 
% dissipation rate and the size of the particles.

function [b_eddy, beta, int_ratio] = BreakageLuoSvendson(iz, im, d, Ns_cell, lambda_min, N_lambdas, params)

    %Debug cases
    if params.break.debug
        t_start.total = cputime;
        fprintf('-Eddy Breakage (Luo): iz = %i; im = %i; d = %0.5f m; m_i = %0.3e;\n', iz, im, d, params.mms(im));
        x = 1;
    end

    %Extract structures
    turb = params.turb;

    %Calculate indices
    xind = (iz-1)*params.Nms + im;

    %Calculate breakage rate based on Luo and Svendson 1996
     if lambda_min(iz) < d %otherwise the bubble cannot be broken by eddies

        %List lambdas and calculate zetas
        lambdas = logspace(log10(lambda_min(iz)), log10(d), N_lambdas);
        zetas = lambdas./d;

        %Constants
        beta_coeff = 2.41;

        %Pull storage vectors
        b_fvd = params.break.bfd_zero;
        ints_fvd = params.break.bfd_zero;

        %Iterate through fvs
        fvhs = (params.fvs_norm(1:end-1) + params.fvs_norm(2:end))./2;
        for iv = 1:length(params.fvs_norm)

            %Pull value
            fv = params.fvs_norm(iv)+1E-12;
            

            %Calculate relevant values
            cf = fv.^(2/3) + (1-fv).^(2/3) -1; %

            %Evaluate integral func
            Xc = @(zeta) -(12.*cf.*params.sigmas(iz))./(beta_coeff .* params.rhos_l(xind) .* turb.eps(iz).^(2/3) .* d.^(5/3) .*zeta.^(11/3));
            integral_func = @(zeta) ((1+zeta).^2)./(zeta.^(11/3)) .* exp(Xc(zeta));
            integral_vals = integral_func(zetas);

            %Calculate breakage rate for this fv
            zeta_int_vals(iv) = trapz(zetas, integral_vals);
            b_fvd(iv) = 0.923 * (1-params.alpha_g(iz)) .* (turb.eps(iz)./d.^2).^(1/3) .* trapz(zetas, integral_vals);

            % %Calculate beta
            % beta(iv) = 2 .* 

        end

        %Calculate overall breakage rate 
        b_eddy = trapz(params.fvs_norm, b_fvd);

        %Calculate bubble size distribution 
        beta = zeta_int_vals./(trapz(params.fvs_norm, zeta_int_vals)); %beta = b_fvd./b_eddy;

        %Expand to full domain
        beta = [beta, fliplr(beta(1:end-1))];

        %Calculate ratio of ints_fvd to beta
        int_ratio = 0;
        % int_ratio = ints_fvd./beta(1:length(ints_fvd));
        % int_ratio = mode(int_ratio);
        % if isnan(int_ratio)
        %     int_ratio = 0;
        % end

        %Check that outputs are physical

     else %The bubble cannot be broken by eddies, so set breakage rate to zero and beta to 1 (uniform distribution)
        b_eddy = 0;
        int_ratio = 0;
        beta = ones(size(params.fvs_norm_all));
     end


end