function [c_src, c_snk, c_mats] = Coalescence(y, params )

    %Make parameters global
    global params

    %Debug message
    if params.debug
        t_start = cputime;
        fprintf('--Coalescence Start\n');
    end

    %Check for inconsistent values
    if any(y < 0)
        warning('Negative numeric densities.')
        y(y < 0) = 0;
    end

    %Allocate vectors for source and sink temrs
    c_src = zeros(length(y), 1); %Consider allocating empty vector in the params to avoid having to make a new one each iteration
    c_snk = c_src;
    c_mats = cell(params.Nz, 1);
    

    %Determien iz inds to use
    if params.sol.sep_layer & strcmp(params.sol.type, 'segregated')
        iz_inds = params.iz;
    else
        iz_inds = 1:params.Nz;
    end
    
    
    %Iterate through all spatial cells
    for iz = iz_inds

        %Pull numeric densities in this spatial celll
    

        cellinds = ((iz-1).*params.Nms + 1):1:(iz.*params.Nms);
        Ns = y();
        Ns_cell = Ns(cellinds);

        %Pull local fluid properties
        uL = 0;


       

        %Calculate average gas velocity from mass dependent velocities
        % V_cell = params.Vcells(iz);
        % V_gas = 0;
        % for ix = 1:params.Nms
        %     V_gas = V_gas + Ns_cell(ix) * V_cell * params.Vms(ix);
        % 
        % end
        % params.alpha_g = V_gas/V_cell;

        %Calculate local fluid properties - only once per cell
        turb = params.turb;

        %Storage matrix for fluxes  between size groups
        c_mat_src = zeros(params.Nms); 
        c_mat_snk = zeros(params.Nms);

        %Storage vector for sink term calculations
        c_ij = zeros(params.Nms);
        cc_ij = c_ij; cw_ij = c_ij; cr_ij = c_ij;



        %Calculate all coalescence rates - needed for both sources/sinks
        for j = 1:params.Nms
            for k = 1:params.Nms
            
                %Pull indices
                %i = xind;
                %j = coalesce_partners(ic,1);
                %k = coalesce_partners(ic,2);
    
                %Caclulate Coalescence Rate
                switch params.coalesce.model
                    case 'Wang_2005'
                        [c_ij(j,k), c_ijs] = Coalescence_Wang(Ns_cell,j,k, iz, params, turb);
                        cc_ij(j,k) = c_ijs(1);
                        cw_ij(j,k) = c_ijs(2);
                        cr_ij(j,k) = c_ijs(3);
                    case 'PrinceBlanch_1990'
                        [c_ij(j,k), c_ijs] = Coalescence_PrinceBlanch(Ns_cell,j,k, iz, params, turb);
                        cc_ij(j,k) = c_ijs(1);
                        cw_ij(j,k) = c_ijs(2);
                        cr_ij(j,k) = c_ijs(3);
                    case 'Scott_1968'
                        c_ij(j,k) = Coalescence_constant(Ns_cell,j,k, iz, params, turb);
                    case 'Hounslow_1988'
                        c_ij(j,k) = Coalescence_constant(Ns_cell,j,k, iz, params, turb);
                end
                
            end
        end

        if any(c_ij(:) < 0)
            x = 1;
        end


        %Iterate through internal variables within the cell
        for ix = 1:(params.Nms-1)

            %Calculate indices
            zind = iz;
            xind = ix; %x refers to bubble size (mass currently)
            yind = (zind-1) * params.Nms + xind;

            %Local density of this mass
            N = y(yind);
            Nb = N .* params.Vsz(zind); %Number of bubbles in the discrete volume

            %Pull coalescence partners and iterate through them
            coalesce_partners = params.coalesce_partners{xind};
            coalesce_bias = params.coalesce_bias{xind};
            coalesce_etas = params.coalesce_etas{xind};
            coalesce_rats = params.coalesce_rats{xind};

            %Handle largest vin
            if ix == params.Nms
               x = 1; 
            end

            %Iterate through 
            cadd = 0;
            csub = 0;
            [N_partners, ~] = size(coalesce_partners);
            for ic = 1:N_partners
        
                %Pull indices
                i = xind; %Index of bin being updated
                j = coalesce_partners(ic,1); %Index of first bubble coalescing
                k = coalesce_partners(ic,2); %Index of second bubble coalescing
                bias = coalesce_bias(ic);
                eta = coalesce_etas(:,ic);
                rat = coalesce_rats(:,ic);
                
                %Don't consider coalescence that creates bubbles larger
                %than upper bound
                if rat > 1
                    rat = 0;
                end
                
                



                %Pull original numeric densities
                jind = (zind-1) * params.Nms + j;
                kind = (zind-1) * params.Nms + k;
                Nj = y(jind);
                Nk = y(kind);

                % %Calculate eta - TEMPORARY
                % if xind < params.Nms
                %     m_jk = params.mms(j) + params.mms(k);
                %     m_low = params.mms(ix-1);
                %     m_mid = params.mms(ix);
                %     m_high = params.mms(ix+1);
                %     if bias < 0
                %         eta_lower = (m_mid - m_jk)./(m_mid - m_low);
                %         eta_higher = (m_jk - m_low)./(m_mid - m_low);
                %     else
                %         eta_lower = (m_high - m_jk)./(m_high - m_mid);
                %         eta_higher = (m_jk - m_mid)./(m_high - m_mid);
                %     end
                % end

                %Caclulate Coalescence Rate - pull from matrix
                c_jk = c_ij(j,k);   %Coalescence_Rate(y,i,j,k, xind, zind, params, turb);

                if c_jk > 0
                    x = 1;
                end

                %Determine rate for each of the three cells
                if bias < 0 %Shares with smaller (left) brackt


                    c_jk_smaller = (1-eta) * c_jk * Nj * Nk;
                    c_jk_middle = eta * c_jk * Nj * Nk;
                    c_jk_larger = 0;

                else %Shares with larger (right) bracket
                    c_jk_smaller = 0;
                    c_jk_middle = eta * c_jk * Nj * Nk;
                    c_jk_larger = (1-eta) * c_jk * Nj * Nk;
                end


                % if bias < 0 %Shares with smaller (left) brackt
                %     c_jk_smaller = eta(1) * c_jk * Nj * Nk;
                %     c_jk_middle = eta(2) * c_jk * Nj * Nk;
                %     c_jk_larger = 0;
                % 
                % else %Shares with larger (right) bracket
                %     c_jk_smaller = 0;
                %     c_jk_middle = eta(1) * c_jk * Nj * Nk;
                %     c_jk_larger = eta(2) * c_jk * Nj * Nk;
                % end

                %Avoid double counting
                if j == k
                    dcoeff = 1; %(1-0.5);
                else
                    dcoeff = 1;
                end
                c_jk_smaller = dcoeff * c_jk_smaller;
                c_jk_middle = dcoeff * c_jk_middle; 
                c_jk_larger = dcoeff * c_jk_larger;

                %Limit rates to be above zero
                c_jk_smaller = max([c_jk_smaller, 0]);
                c_jk_middle = max([c_jk_middle, 0]);
                c_jk_larger = max([c_jk_larger, 0]);
                
                
                %Store results
                c_src(yind) = c_src(yind) + c_jk_middle;
                if ix > 1
                    c_src(yind-1) = c_src(yind-1) + c_jk_smaller;
                end
                if ix < params.Nms
                    c_src(yind+1) = c_src(yind+1) + c_jk_larger;
                end

                %Store fluxes between representative sizes
                i_smaller = i-1; i_smaller = max([i_smaller,0]);
                i_middle = i;
                i_larger = i+1; i_larger = min([i_larger, params.Nms]);
                mjk = params.mms(j) + params.mms(k);
                Xj = params.mms(j)/mjk; %Can be calculated beforehand to reduce cost (marginally)
                Xk = params.mms(k)/mjk;
                c_mat_src(j, i_smaller) = c_mat_src(j, i_smaller) + c_jk_smaller * Xj;%Mass frrom smaller bubble
                c_mat_src(j, i_middle) = c_mat_src(j, i_middle) + c_jk_middle * Xj;
                c_mat_src(j, i_larger) = c_mat_src(j, i_larger) + c_jk_larger * Xj;%Mass from larger bubble
                c_mat_src(k, i_smaller) = c_mat_src(k, i_smaller) + c_jk_smaller * Xk;%Mass frrom smaller bubble
                c_mat_src(k, i_middle) = c_mat_src(k, i_middle) + c_jk_middle * Xk;
                c_mat_src(k, i_larger) = c_mat_src(k, i_larger) + c_jk_larger * Xk;%Mass from larger bubble
             
                %Subtract bubbles from their original bins
                if j == k
                    c_snk(jind) = c_snk(jind) + 2*c_jk * Nj * Nk; 
                else
                    c_snk(jind) = c_snk(jind) + c_jk * Nj * Nk;
                    c_snk(kind) = c_snk(kind) + c_jk * Nj * Nk; 
                end

                %
                %c_mat_snk(j) = c_snk(jind)
                %c_mat_snk(j)



                %Calculate total change for this bracket
               % c_src(i) = c_src(i) + c_jk;

                %Log result for sink calculations
               % c_ij(j,k) = c_jk;
    
            end

            

        end

        %c_snk(i) = c_snk(i) - 0; %sum(Ns .* );


        if any(c_mat_src(:) > 0)
            x = 1;
        end

        %Store c_mat
        c_mat.src = c_mat_src;
        c_mat.snk = c_mat_snk;
        c_mats{iz} = c_mat_src;


    end

    %Catch negative values
    if any(c_src < 0)
        x = 1;
    end

    %Check mass conservation and rescale 
    if params.src.debug
        N_src = sum(c_src);
        N_snk = sum(c_snk);
        params.coalesce.m_src(params.src.its) = sum(c_src .* [repmat(params.mms, 1, params.Nz)]'); %params.mms_rep); %[repmat(params.mms, 1, params.Nz)]');
        params.coalesce.m_snk(params.src.its) = sum(c_snk .* [repmat(params.mms, 1, params.Nz)]'); %params.mms_rep); %[repmat(params.mms, 1, params.Nz)]');
        c_src = (params.coalesce.m_snk(params.src.its)./(params.coalesce.m_src(params.src.its) + 1E-16)) * c_src;
    
    end


    
    %m_src = sum(c_src .* [repmat(params.mms, 1, params.Nz)]');
    %m_snk = sum(c_snk .* [repmat(params.mms, 1, params.Nz)]');

    % if sum(m_src) > 0
    %     c_src = (sum(m_snk)./sum(m_src)) .*c_src;
    % end

    %Debug message
    if params.debug
        t_end = cputime;
        t_req = t_end - t_start;
        fprintf('--Coalescence End (t = %0.4f s)\n', t_req);
    end



end





function [c_jk, c_jks] = Coalescence_Wang(y, j,k, zind, params, turb)

    %Define index to sample from for liquid properties
    i = (zind - 1) * params.Nms + 1;

    %Pull mass values at pivots
        %mi = params.mms(i);
    mj = params.mms(j);
    mk = params.mms(k);

    %Determine which bracket 

    %Pull other values at pivots
    jind = j; %jind = (zind - 1) .* params.Nms + j;
    kind = k; %(zind - 1) .* params.Nms + k;
    Nj = y(jind);
    Nk = y(kind);
    %Dd_j = ((6 .* params.mbs(j+1) .* params.R .* params.Ts(i))./(pi .* params.ps(i) .* params.Ms(i))).^(1/3)...
      %  - ((6 .* params.mbs(j) .* params.R .* params.Ts(i))./(pi .* params.ps(i) .* params.Ms(i))).^(1/3);
    %Dd_k = ((6 .* params.mbs(k+1) .* params.R .* params.Ts(i))./(pi .* params.ps(i) .* params.Ms(i))).^(1/3)... 
      %  - ((6 .* params.mbs(k) .* params.R .* params.Ts(i))./(pi .* params.ps(i) .* params.Ms(i))).^(1/3);

                % %Pull mass values of neighboring brackets
                % mi_lo = params.mms(i-1);
                % if i == params.Nms
                %     mi_hi = params.mbs(i+1);
                % else
                %     mi_hi = params.mms(i+1);
                % end

    %Calculate coalesced mass
    mjk = mj + mk;

    %Calculate diameters of bubbles 
        %di = ((6 * mi * params.R .* params.Ts(i))./(params.ps(i) .* params.Ms(i))).^(1/3);
    dj = params.d_mu(zind, jind); %((6 * mj * params.R .* params.Ts(i))./(pi .* params.ps(i) .* params.Ms(i))).^(1/3);
    dk = params.d_mu(zind, kind); %((6 * mk * params.R .* params.Ts(i))./(pi .* params.ps(i) .* params.Ms(i))).^(1/3);

    %Coalescence due to turbulent eddies - Adapted from Wang et al. 2005 -
    %   Theoreetical prediction of flow regime ...
    u_bar_j = sqrt(2) .* (turb.eps(zind) .* dj).^(1/3); %Mean turbulent velocity 
    u_bar_k = sqrt(2) .* (turb.eps(zind) .* dk).^(1/3);
    nj = Nj/params.dbds(zind, j); %params.mds(j);%dbds(j); %CHECK CHECK CHECK
    nk = Nk/params.dbds(zind, k); % params.mds(k); %dbds(k);  
    hb_jk = 6.3 .* (nj + nk).^(-1/3); 
    lbt_j = sqrt(2) * (turb.eps(zind).*dj).^(1/3) .* (((dj/2).^2)/turb.eps(zind)).^(1/3);
    lbt_k = sqrt(2) * (turb.eps(zind).*dk).^(1/3) .* (((dk/2).^2)/turb.eps(zind)).^(1/3);   
    lbt_jk = sqrt(lbt_j^2 + lbt_k^2); %lbt_jk = sqrt( (0.89 .* dj).^2 + (0.89 .* dk).^2);
    gamma_jk = exp(-(hb_jk./lbt_jk).^6);
    omega_c = (pi./4) .* ((params.alpha_g_max)./(params.alpha_g_max - params.alpha_g(zind))) .* gamma_jk...
        .* sqrt(2) .* turb.eps(zind).^(1/3) .* (dj + dk).^2 .* sqrt(dj.^(2/3) + dk.^(2/3));  %sqrt(u_bar_j.^2 + u_bar_k.^2);
    

    %omega_c = (pi/4) .* (dj + dk).^2 .* sqrt(u_bar_k.^2 + u_bar_j.^2);


    psi = 1;
    u_bar_jk = sqrt(u_bar_j.^2 + u_bar_k.^2);
    We_ij = (params.rhos_l(i) .* (dj) .* u_bar_jk.^2)./params.sigmas(zind);
    rho_g = 0.5;
    gamma_VM = 0.5; %coefficient of virtual mass
    Pc_jk = exp(-psi .* sqrt(0.75 .* (1+ (dj/dk).^2) .* (1 + (dj./dk).^3))./((rho_g/params.rhos_l(i) + gamma_VM) .* (1 + (dj./dk)).^3) .* sqrt(We_ij));
    cc_jk = omega_c * Pc_jk;

    if gamma_jk > 0
        x = 1;

    end


    if params.coalesce.rise
        %
        uj = params.uz_mu(zind, jind);
        uk = params.uz_mu(zind, kind);
    
        Pr_jk = 0.5; %Efficiency of coalescence for j and k due to different rise velocities
        cr_jk = (pi/4) .* ((dj + dk).^2) .* abs(uj - uk)*Pr_jk;
    else
        cr_jk = 0;
    end


    %Coalescence due to wake entrainment
    K_w1 = 15.4;
    K_w2 = 0.46;
    dc_w = 4 .* sqrt(params.sigmas(zind)./(params.g .* params.rhos_l(i))); %Critica diameter for wake entrainment
    if dj > (dc_w/2) || dk > (dc_w/2)
        
        if dj  < dk
    
            %Calculate slip velocity - seems to be for a cap
            %bubble
            u_bar_slip = 0.71 .* sqrt(params.g.*dk);
    
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
            cw_jk = K_w1 .* theta_c .* dk.^2 .* u_bar_slip .* exp(-K_w2 .* (((params.rhos_l(zind)^3 .* turb.eps(zind).^2)./(params.sigmas(zind).^3) .* (dj .* dk./(dj + dk)).^5).^(1/6)));    %(params.rhos_l(i).^(0.5) * turb.eps(zind).^(1/3))/sqrt(params.sigmas(i)) .* (dj.*dk./(dj +dk)).^(5/6)); %exp(-K_w2 .* (((params.rhos_l(i).^(3) .* turb.eps(zind).^2)./(params.sigmas(zind).^3)) .* ((dj .* dk)./(dj + dk)).^5).^(1/6));
    
        else % j is the larger or equivalent bubble
    
            %Calculate slip velocity
            u_bar_slip = 0.71 .* sqrt(params.g.*dj);
    
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
            cw_jk = K_w1 .* theta_c .* dj.^2 .* u_bar_slip .* exp(-K_w2 .* (((params.rhos_l(zind)^3 .* turb.eps(zind).^2)./(params.sigmas(zind).^3) .* (dj .* dk./(dj + dk)).^5).^(1/6)));    %(params.rhos_l(i).^(0.5) * turb.eps(zind).^(1/3))/sqrt(params.sigmas(i)) .* (dj.*dk./(dj +dk)).^(5/6)); %exp(-K_w2 .* (((params.rhos_l(i).^(3) .* turb.eps(zind).^2)./(params.sigmas(zind).^3)) .* ((dj .* dk)./(dj + dk)).^5).^(1/6));
    %cw_jk = K_w1 .* theta_c .* dj.^2 .* u_bar_slip .* exp(-K_w2 .* (params.rhos_l(i).^(0.5) * turb.eps(zind).^(1/3))/sqrt(params.sigmas(i)) .* (dj.*dk./(dj +dk)).^(5/6)); %exp(-K_w2 .* (((params.rhos_l(i).^(3) .* turb.eps(zind).^2)./(params.sigmas(zind).^3)) .* ((dk .* dj)./(dk + dj)).^5).^(1/6));
        end
    else
        cw_jk = 0;
    end

    %Scale wake down
    %cw_jk = 0.006 * cw_jk; %TEMPORARY



    %Calculate probability for bubbles to coalesce 
    c_jks = [cc_jk, cw_jk, cr_jk];

    %Calculate coalescence rate total
    cjk = 0;
    if params.coalesce.eddy
        cjk = cjk + cc_jk;
    end
    if params.coalesce.wake
        cjk = cjk + cw_jk;
    end
    if params.coalesce.rise
        cjk = cjk + cr_jk;
    end
    
    %cjk = cc_jk + cw_jk; % Coalescence rate #/m^3-s - Used direclty for coalescence sink term



    %Calculate distribution coeff - MOVE OUTSIDE INTEGRATION
    %LOOP - PUT IN PBM
    % if mjk < mi
    %     eta_ijk = (mjk - mi_lo)./(mi - mi_lo);
    % elseif mjk >= mi
    %     eta_ijk = (mi_hi - mjk)./(mi_hi - mi);
    % end

    %Calculate total change for this bracket
    %c_src(i) = c_src(i) + dcoeff * eta_ijk * cjk * y(jind) * y(kind);
    

    c_jk = cjk;


    %c_jk = dcoeff * eta_ijk * cjk * y(jind) * y(kind); %USE FOR
    %DISTRIBUTIN

      %  (1 - 0.5 * )
    %eta_ijk = 1;

    if c_jk < 0
        x = 1;
    end



end


function cc_jk = Coalescence_Turb(eps_turb, dj, dk)

    u_bar_j = sqrt(2) .* (eps_turb .* dj).^(1/3);
    u_bar_k = sqrt(2) .* (eps_turb .* dk).^(1/3);
    hb_jk = 6.3 .* (Nj/Dd_j + Nk/Dd_k).^(-1/3);
    lbt_jk = sqrt( (0.89 .* dj).^2 + (0.89 .* dk).^2);
    gamma_jk = exp(-(hb_jk./lbt_jk).^6);
    omega_c = (pi./4) .* ((alpha_g_max)./(alpha_g_max - alpha_g)) .* gamma_jk...
        .* sqrt(2) .* eps_turb.^(1/3) .* (dj + dk).^2 .* sqrt(dj.^(2/3) + dk.^(2/3));  %sqrt(u_bar_j.^2 + u_bar_k.^2);
    psi = 1;
    u_bar_jk = sqrt(u_bar_j.^2 + u_bar_k.^2);
    We_ij = (params.rhos_l(i) .* ((dj + dk)./2) .* u_bar_jk.^2)./params.sigmas(i);
    rho_g = 0.5;
    Pc_jk = exp(-psi .* sqrt(0.75 .* (1+ (dj/dk).^2) .* (1 + (dk./dk).^3))./((rho_g/params.rhos_l(i) + 0.5) .* (1 + (dj./dk)).^3) .* sqrt(We_ij));
    cc_jk = omega_c * Pc_jk;


end


function Coalescence_Wake



end





function c_jk = Coalescence_constant(y, j,k, zind, params, turb)

    c_constant = params.coalesce.constant_rate; %Applies to all 
    

    c_jk = 0.5*c_constant;

end



function c_jk = Coalescence_Lehr(y, j,k, zind, params, turb)

    %Define index to sample from for liquid properties
    i = (zind - 1) * params.Nms + 1;

    %Pull mass values at pivots
    mj = params.mms(j);
    mk = params.mms(k);

    %Pull other values at pivots
    jind = j; %jind = (zind - 1) .* params.Nms + j;
    kind = k; %(zind - 1) .* params.Nms + k;
    Nj = y(jind);
    Nk = y(kind);

    %Calculate coalesced mass
    mjk = mj + mk;

    %Calculate diameters of bubbles 
    dj = params.d_mu(zind, jind); 
    dk = params.d_mu(zind, kind);

    




end