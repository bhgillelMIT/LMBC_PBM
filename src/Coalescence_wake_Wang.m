function cw_jk = Coalescence_wake_Wang(dj, dk, turb, zind, i, params)

    %Coalescence due to wake entrainment
    K_w1 = 1.2; %0.151; %15.4;
    K_w2 = 0.46;
    dc_w = 4 .* sqrt(params.sigmas(zind)./(params.g .* params.rhos_l(i))); %Critica diameter for wake entrainment
    L_w_d_ratio = 3;
    kappa = 1;
    k1 = 1;
    k_empirical = 1;
    if dj > (dc_w/2) || dk > (dc_w/2)
        
        if dj  < dk
    
            %Calculate coefficient in detail
            L_w = L_w_d_ratio * dk;
            K = params.coalesce.c_wake; %(pi/4) * k_empirical./((L_w/dk) - 0.5) .* ((L_w/(dk/2))^(1/3) - 1);

            %Calculate slip velocity - seems to be for a cap
            %bubble
            u_bar_slip = params.uz_mu(params.xinds(i)); %u_bar_slip = 0.71 .* sqrt(params.g.*dk);
    
            %Determine if bubble is above critical size for wake
            %entrainment
            if dk >= dc_w/2
                theta_c = ((dk - dc_w/2).^6)./((dk - dc_w/2).^6 + (dc_w/2).^6);
            else
                theta_c = 0;
            end
    
            %Calculate coalescence rate
           % omega_w = K .* dk.^2 .* u_bar_slip;
            %cw_jk = 0.0073 * u_bar_slip * theta_c * (params.alpha_g(zind).^2)./dk.^4;
            Pw_jk = exp(-K_w2 .* (((params.rhos_l(zind)^3 .* turb.eps(zind).^2)./(params.sigmas(zind).^3) .* (dj .* dk./(dj + dk)).^5).^(1/6)));
            cw_jk = K .* theta_c .* dk.^2 .* u_bar_slip .* Pw_jk;    %(params.rhos_l(i).^(0.5) * turb.eps(zind).^(1/3))/sqrt(params.sigmas(i)) .* (dj.*dk./(dj +dk)).^(5/6)); %exp(-K_w2 .* (((params.rhos_l(i).^(3) .* turb.eps(zind).^2)./(params.sigmas(zind).^3)) .* ((dj .* dk)./(dj + dk)).^5).^(1/6));
    
        else % j is the larger or equivalent bubble

            %Calculate coefficient in detail
            L_w = L_w_d_ratio * dj;
            K = params.coalesce.c_wake; %(pi/4) * k_empirical./((L_w/dj) - 0.5) .* ((L_w/(dj/2))^(1/3) - 1);
    
            %Calculate slip velocity
            u_bar_slip = params.uz_mu(params.xinds(i)); %u_bar_slip = 0.71 .* sqrt(params.g.*dj);
            u_bar_slip_alt = 0.71 .* sqrt(params.g.*dj);
    
            %Determine if bubble is above critical size for wake
            %entrainment
            if dj >= dc_w/2
                theta_c = ((dj - dc_w/2).^6)./((dj - dc_w/2).^6 + (dc_w/2).^6);
            else
                theta_c = 0;
            end
    
            %Calculate coalescence rate
           % omega_w = K_ .* dj.^2 .* u_bar_slip;
            %cw_jk = 0.0073 * u_bar_slip * theta_c * (params.alpha_g(zind).^2)./dj.^4;
            cw_jk = K .* theta_c .* dj.^2 .* u_bar_slip .* exp(-K_w2 .* (((params.rhos_l(zind)^3 .* turb.eps(zind).^2)./(params.sigmas(zind).^3) .* (dj .* dk./(dj + dk)).^5).^(1/6)));    %(params.rhos_l(i).^(0.5) * turb.eps(zind).^(1/3))/sqrt(params.sigmas(i)) .* (dj.*dk./(dj +dk)).^(5/6)); %exp(-K_w2 .* (((params.rhos_l(i).^(3) .* turb.eps(zind).^2)./(params.sigmas(zind).^3)) .* ((dj .* dk)./(dj + dk)).^5).^(1/6));
    %cw_jk = K_w1 .* theta_c .* dj.^2 .* u_bar_slip .* exp(-K_w2 .* (params.rhos_l(i).^(0.5) * turb.eps(zind).^(1/3))/sqrt(params.sigmas(i)) .* (dj.*dk./(dj +dk)).^(5/6)); %exp(-K_w2 .* (((params.rhos_l(i).^(3) .* turb.eps(zind).^2)./(params.sigmas(zind).^3)) .* ((dk .* dj)./(dk + dj)).^5).^(1/6));
        end
    else
        cw_jk = 0;
    end




end