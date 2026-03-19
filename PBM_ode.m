function dydt = PBM_ode(t,y, params)

    %Make parameters global
    global params

    %Run parent function if called directly
    if nargin == 0
        PBM_v3();
    end

    %Settings
    debug = params.debug;

    %y = real(y);

    %Model constants
    %alpha_g_max = 0.8; %Maximum gas holdup


    %Pull parameters
    %uz = params.uz;
    bottom_inds = params.bottom_inds;
    top_inds = params.top_inds;

    %Normalize numeric densities if 
    if params.sol.single_layer
        m_total = sum(params.mms_rep .* y .* params.Vcells_rep);
        m_ref = params.m_total(1);
        y = y .* (m_ref/m_total);
    end

    %Isolate variables
    dFdt = zeros(params.N_volumes, 1);
    F = y; %F(F < 0) = 0;
    Fs = F; Fs(F < 0) = 0;

    %Check for nans
    if any(isnan(y))
        warning('NaNs');
    end


    %Calculate mean temperature and conversion for each representative mass
    params = CalcTempConvMeans(Fs, params);

    %Calculate total numeric density for each 
    [params.Ns_z, params.Ns_m, params.Ns_T, params.Ns_fracs] = CalcNumericDensities(y, params);

    %Identify indices for inlet
    inlet_inds = find(params.zinds == 1); 
    
    %Warnings
    if any(isnan(F))
        warning('NaNs');
    end
    if any(y < 0)
        warning('Negatives');
    end

    %Calculate local conditions
    Xs = params.Xs;
    Ts = params.Ts;
    ns_CH4 = params.ns_i .* (1 - Xs); 
    ns_H2 = 2 * (params.ns_i - ns_CH4);
    Ms = (ns_CH4.*params.M.CH4 + ns_H2.*params.M.H2)./(ns_CH4 + ns_H2); %params.Xs
    Vs = (params.ms .* params.R .* Ts)./(params.ps .* Ms);
    Ds = (6 .* Vs./pi).^(1/3);
    %Eobs = (params.tin.density(Ts) .* params.g .*Ds.^2 )./params.tin.surf_tension(Ts);
    %uzs = params.ub_func(Ts, Ds); %CalcVelocity(Ds, Eos, params);

    %Log results
    %params.uzs = uzs;
    params.Ms = Ms;

    %Calculate local hydrodynamic conditions - location dependent
    

    %Calculate flux limiters
    %phis = vanleer(F)

    % Compute slopes using MUSCL with Van Leer limiter
    rz = zeros(1, length(F));
    rx = zeros(1, length(F));
    for j = 1:length(params.Ns_m)

        %Pull corresponding indices
        xind = mod(j-1, params.Nms)+1;
        zind = floor((j-1)/params.Nms) + 1;

        %Determine adjacent indices
        ind_above = j - params.Nms;
        ind_below = j + params.Nms;

        %Calcualte rz
        if zind == 1 | zind == params.Nz %Bottom or top layer
            rz_val = 0;
        else
            rz_val = (F(j) - F(ind_below))/(F(ind_above) - F(j) + 1E-12);
        end

        %Determine which indices to apply this to
        relinds = find(params.zinds == zind & params.xinds == xind);

        %Apply value to 
        rz(relinds) = rz_val;

        
        
        %Identify the index of the matching cell upstream
        % if zind > 1
        % 
        %     %Identify upstream T index
        %     merged = params.T_z.chars.merged{xind};
        %     T_match = merged.bins_in(:,zind) == Tind; %Need to handle multiple incoming cells
        % 
        % 
        %     %Identify other upstream indices
        %     zmatch = params.zinds == (zind-1);
        %     xmatch = params.xinds == xind;
        % 
        %     %Identify incoming 
        %     Tmatch = 1;%Complicated by merging 
        %     Xmatch = params.Xinds == Xind;
        %     upstream_ind = find(zmatch & xmatch & Tmatch & Xmatch);
        %     downstream_ind = 1;
        % 
        %     %NOTE: Only consider the total numeric density under each mass
        %     %bin for calculating the limiter 
        % 
        % 
        % end

            
        

        %Calculate rx
        % if xind == 1 | xind == params.Nms
        %     rx(j) = 0;
        % else
        %     rx(j) = (F(j) - F(j-1))/(F(j+1) - F(j) + 1E-12);
        % end

        
        % %Determine indices
        % xind = params.xinds(j);
        % Tind = params.Tinds(j);
        % Xind = params.Xinds(j);
        % cellind = params.cellinds(j);
        % zind = params.zinds(j);


    end

    %Apply Van Leer Limiter
    %phix = (rx + abs(rx)) ./ (1 + abs(rx)); phix = phix'; % Van Leer limiter
    phiz = (rz + abs(rz)) ./ (1 + abs(rz)); phiz = phiz';

    %Calculate derivatives using MUSCL scheme


    % %Calcualte edge values using MUSCL scheme
    % if strcmp(params.scheme, 'MUSCL')
    %     FR_left = F(1:end-2) + 0.5 .* phi(1:end-2) .* (F(2:end-1) - F(1:end-2));
    %     FR_right = F(2:end-1) - 0.5 .* phi(2:end-1) .* (F(3:end) - F(2:end-1));
    %     FR = (FR_left + FR_right)./2; FR = [FR;(F(end) + F(end-1))./2;F(end)];
    % 
    %     FL_left = F(1:end-2) + 0.5 .* phi(1:end-2) .* (F(2:end-1) -  F(1:end-2));
    %     FL_right = F(2:end-1) - 0.5 .* phi(2:end-1) .* (F(3:end) - F(2:end-1));
    %     FL = (FL_left + FL_right)./2; 
    %     FL = [0; FL; (F(end-1) + F(end))./2];
    % else
    %     phi = ones(size(F));
    %     FL = F(1:end-1) + 0.5 .* phi(1:end-1) .* (F(2:end) - F(1:end-1)); % Left-biased
    %     FR = F(2:end) - 0.5 .* phi(2:end) .* (F(2:end) - F(1:end-1));
    %     FL = [0;FL]; FR = [FR;0];
    % end



    %Calculate current gas holdup in each spatial cell
    params = CalcLocalProperties(Fs, params); %gas fraction, turbulence, etc.   %params.uzs = 1;

    %Calculate source terms
    if t < params.sol.src_delay
        h = zeros(size(F));
    else
        params.src.its = params.src.its + 1;
        cadd = zeros(length(F), 1); csub = cadd;
        badd = cadd; bsub = cadd;
        if params.coalesce.active
            [cadd, csub, cmats] = Coalescence(params.Ns_m, params); %Coalescence(Fs, params );  
            cadd = cadd * params.coalesce.damper; %Damper is only used for debugging
            csub = csub * params.coalesce.damper;
            params.cmats = cmats; %Store for heat/reaction calculations
        end
        if params.break.active
            if params.sol.solve_details
                [badd, bsub, bmats] = Breakage(params.Ns_m, params); %Breakage(Fs, params);
                params.break.badd = badd;
                params.break.bsub = bsub;
                params.bmats = bmats;
            else
                badd = params.break.badd;
                bsub = params.break.bsub;
            end
            badd = params.break.damper * badd;
            bsub = params.break.damper * bsub;
        end

        %Distribute
        if params.heat.active
            h = params.h;

            h_m = cadd(:) - csub(:) + badd(:) - bsub(:);

            %Store in params
            params.cadd = cadd; params.csub = csub;
            params.badd = badd; params.bsub = bsub;

            %Function to distribute proportionally
            h = DistSourceTerms(y, h_m, params);


            %Iterate through bins and allocate
            ind = 1; %index for coalescence and breakage vectors 
            for iz = 1:params.Nz
                for im = 1:params.Nms
                    subinds = find(params.zinds == iz & params.xinds == im);
                    h(subinds) = params.Ns_fracs(subinds) .* (cadd(ind) - csub(ind)...
                        + badd(ind) - bsub(ind));
                    ind = ind + 1;
                end
            end
        else
            h = cadd(:) - csub(:) + badd(:) - bsub(:);
        end

         
    end

    if params.sol.single_layer
        h = h;
    end

    if h(1) == 10000
        figure();
        subplot(1,2,1);
        plot(y);
        subplot(1,2,2);
        plot(cadd); hold on;
        plot(csub); plot(badd); plot(bsub);
        legend('cadd', 'csub', 'badd', 'bsub');
        x = 1;
    end

    % %Fade source term in
    % if t < params.sol.src_delay
    %     h = h .* t/params.sol.src_delay;
    % end

    %Constrain source terms to avoid negative numeric densities
    delta_max = h * params.dt_max;
    F_min = F + delta_max;
    F_ratio = zeros(size(F));
    neg_inds = F_min < 0;
    F_ratio(neg_inds) = delta_max(neg_inds)./F(neg_inds);
    %F_ratio = delta_max./F; 
    %F_ratio(isinf(F_ratio)) = 0;

    

    % if any(F_min < 0)
    %     neg_inds = F_ratio < 0;
    %     h_adj_fac = 0.98 .* double(1)./double(F_ratio(neg_inds));
    %     if ~isempty(h_adj_fac)
    %         h(neg_inds) = h_adj_fac .* h(neg_inds);
    %     end
    % end
    % 
    % if any(isnan(h))
    %     warning('NaNs');
    % end

    if params.sol.single_layer

        %Assign changes
        dFdt = h;

        %Print update
        if debug
            fprintf('PBM (t = %0.6f s)\n', t)
        end

    else

        %Iterate through cells and impose equations
        FLs = zeros(size(F));
        FRs = zeros(size(F));
        for i = 1:length(y)
    
            %Pull cell
            zind = params.zinds(i);
            xind = params.xinds(i);
            Tind = params.Tinds(i);
            Xind = params.Xinds(i);
            cellind = params.cellinds(i);
            volcell = params.mesh.volcells{cellind}; %spatial cell 
    
            %Determine cell 
    
            %Calculate derivative
            %stencil_vals = params.FD_mat(i,:);
            %stencil_vals = stencil_vals(stencil_vals ~= 0);
            stencil_vals = volcell.cent_FD_coeffs_D1(:,2);
            stencil_inds = (volcell.stencil_cells(:,2) - 1) * params.Nms + xind;
            dFdz = F(stencil_inds)' * stencil_vals;
    
            %MUSCL for z
            %phiz(i) = 1;
            switch params.sol.scheme
                case 'MUSCL'
    
                    if zind == 1
                        zi_up = i + params.Nms;
                        zi_upup = i + 2*params.Nms;
                        zi_down = i; %No z below
                        FR_left = F(i) + 0.5 * phiz(i) * (F(zi_up) - F(i));
                        FR_right = F(zi_up) - 0.5 * phiz(zi_up) * (F(zi_upup) - F(zi_up));
                        FR = FR_left; %(FR_left + FR_right)/2;
                        FL = FR; %Boundary condition - dirichlet
                    elseif zind == params.Nz
                        zi_down = i - params.Nms;
                        FL = F(zi_down); % + 0.5*phiz(zi_down)*(F(i) - F(zi_down));
                        FR = F(i);
                    elseif zind == params.Nz - 1
                        zi_up = i + params.Nms;
                        zi_down = i - params.Nms;
                        
                        FL_left = F(zi_down) + 0.5*phiz(zi_down)*(F(i) - F(zi_down));
                        FL_right = F(i) - 0.5*phiz(i)*(F(zi_up) - F(i));
                        FL = FL_left; %(FL_left + FL_right)/2;
                        
                        FR = F(i) + 0.5 * phiz(i) * (F(zi_up) - F(i));
                    else
                        zi_up = i + params.Nms;
                        zi_upup = i + 2*params.Nms;
                        zi_down = i - params.Nms;
                        FR_left = F(i) + 0.5 * phiz(i) * (F(zi_up) - F(i));
                        FR_right = F(zi_up) - 0.5 * phiz(zi_up) * (F(zi_upup) - F(zi_up));
                        FR = (FR_left + FR_right)/2;
        
                        FL_left = F(zi_down) + 0.5*phiz(zi_down)*(F(i) - F(zi_down));
                        FL_right = F(i) - 0.5*phiz(i)*(F(zi_up) - F(i));
                        FL = FL_left; %(FL_left + FL_right)/2;
                    end
        
                    %Log results
                    FLs(i) = FL;
                    FRs(i) = FR;
    
                case 'Upwind'
                    if zind == 1
                        FR = F(i);
                        FL = FR;
                    elseif zind == params.Nz
                        FL = sum(F(params.inds_below{i})); %F(i-params.Nms);
                        FR = F(i);
                    else
                        FL = sum(F(params.inds_below{i})); %F(i-params.Nms);
                        FR = F(i);  
                    end
                case 'None'
                    FR = 0; FL = 0;
    
            end
    
            %Debug
            if zind > 1 & FL > 0 & FR > 0 & t > 2.5
               x = 1;
            end
    
            %Calculate local properties
            
    
    
    
            %Calculate x derivative
            xc = params.rmesh.xsc(xind);
            x_stencil = params.rmesh.stencil_inds_c(:,xind) + (cellind - 1) * params.Nr;
            x_coeffs = params.rmesh.coeffs_c(:,xind);
            dFdx = F(x_stencil)' * x_coeffs;
    
            %Calculate other derivatives
            dpdt = -params.reactor.rho_L_bar * params.g * params.ug(zind);
            dxdp = -(1/3) .* ((0.75 * params.nRTs(xind))/(pi * params.ps(i)^4)).^(1/3);
            dxdt = dxdp * dpdt;
            %dxdt = -dpdt .* (params.nRTs(xind)./(36 * pi * params.ps(i))).^(1/3);
            d2xdxdt = (params.reactor.rho_L_bar * params.g * params.ug(zind))/3 ...
                        * (params.p_orifice/(params.ps(i)^4))^(1/3); %params.rho_L * params.g * params.uz * ((params.p_orifice)./(9*params.ps(i)))^(1/3);
    
            %Calculate coefficients
            uz = params.uz_mu(zind, xind);
            dmdt = 0; 
            cx = dxdt;
            cz = uz;
            cF = d2xdxdt;
    
            %Calculate more accurate velocity estimate at inferfaces
    
    
            
            %Debug
            if t > 0.01
                x = 1;
            end
    
            %Handle boundary condition
            if any(cellind == bottom_inds)
                if params.sol.orifice_BC_type == 1
                    dFdt_in = zeros(1,params.Nms);
                elseif params.sol.orifice_BC_type == 2
                    dFdt_in = (params.N_dot_os(i) - (uz * FR))/params.dz;
                end
            end
    
            %Solve equation
            switch params.sol.scheme
                case 'MUSCL'
                    stencil_vals = volcell.cent_FD_coeffs_D2(:,2);
                    stencil_inds = (volcell.stencil_cells(:,2) - 1) * params.Nms + xind;
                    d2Fdz2 = F(stencil_inds)' * stencil_vals;
                    if any(cellind == bottom_inds)
                        dFdt(i) = dFdt_in; %((uz * FL) - (uz * FR))/params.dz + params.mu_art * d2Fdz2;
                    else
                        dFdt(i) = h(i) + ((uz * FL) - (uz * FR))/params.dz + params.sol.mu_art * d2Fdz2;
                    end
                case 'Upwind'
                    if any(cellind == bottom_inds)
                        dFdt(i) = dFdt_in; %((uz * FL) - (uz * FR))/params.dz;
                    else
                        dFdt(i) = h(i) + ((uz * FL) - (uz * FR))/params.dz;
                    end
                case 'None'
                    dFdt(i) = h(i);
           
                otherwise
                    if any(cellind == bottom_inds)
                        dFdt(i) = 0;
                    %elseif t < 0.01
                    %    dFdt(i) = h - cz *dFdz;
                    else
                        %dFdt(i) = h -cz * dFdz - cx * dFdx - cF * F(i);
                        dFdt(i) = h(i) -cz * dFdz;
                    end
            end
    
            %Debug
            if dFdt(i) > 0
                x=1;
            end
    
            
    
    
        end
        
        %Create debug table
        if debug
            fprintf('PBM (t = %0.6f s)\n', t)
            % debugtable = table(F, h, dFdt);
    
            % figure();
            % plot(y); hold on;
            % plot(dFdt);
            % plot(h);
            % plot(cadd);
            % plot(badd);
            % plot(-csub);
            % plot(-bsub);
            % 
            % 
            % close all
        end
    end


    %Limit dFdt
    % dF_max = dFdt * params.dt_max;
    % F_min = F + dF_max;
    % neg_inds = F_min < 0;
    % F_ratio = zeros(size(F));
    % nonzeros = F ~= 0;
    % F_ratio(nonzeros) = abs(dF_max(nonzeros)./F(nonzeros)); 

    
    
    % 
    % if any(F_min < 0)
    %     dFdt(F_min < 0) = -0.98 * F(F_min < 0)/params.dt_max; %
    %     x = 1;
    % end



    % %
    % dFdt(F < 0 & dFdt < 0) = 0;

    %Define output vectro
    dydt = dFdt;



end


function phi = vanleer(F)

    %Calculate gradient ratio
    r = zeros(size(F));
    for i = 2:length(F) - 1
        if (F(i+1) - F(i)) ~= 0
            r(i) = (F(i) - F(i-1))./(F(i+1) - F(i));
        else
            r(i) = 0;
        end
    end

    %Calculate phi
    phi = (r + abs(r))./(1 + abs(r));


end


function flux(F, x)


end


function dFdx = Calc_dFdx()

    

end



function ub = CalcVelocity(Db, rho_L, mu_L, sigma, Eob)

    %Physical constants
    g = 9.81;

    %Initial velocity guess
    ub_guess = 0.3;

    ub_func = @(ub) sqrt((4*Db*g)./(3 * CalcCd(ub, Db, rho_L, mu_L, sigma, Eob))) - ub;
    ub = fsolve(ub_func, ub_guess);
    
    


end

function Cd = CalcCd(ub, Db, rho_L, mu_L, sigma, Eob)
    Reb = (rho_L * ub * Db)/mu_L;
    Cd_sph = (24/Reb)*(1 + 0.15*Reb^(0.687));
    Cd_ell_cap = 8/3 * Eob/(Eob + 4);
    Cd = max([Cd_sph, Cd_ell_cap]);

end


% %Calculate RHS of PBM
% 
% 
% 
% function params = CalcLocalProperties(y, params)
% 
%     %Resolve local gas velocity
%     %params.ug = 0.3; %TEMPORARY
% 
%     %Iterate through spatial cells
%     for iz = 1:params.Nz
% 
%         %Pull current z
%         z = params.zms(iz);
% 
%         %Pull numeric densities in this spatial celll
%         cellinds = ((iz-1).*params.Nms + 1):1:(iz.*params.Nms);
%         Ns_cell = y(cellinds);
% 
% 
%         %Calculate gas holdup
%         V_cell = params.Vcells(iz);
%         V_gas = 0;
%         V_dot = zeros(1, params.Nms);
%         for ix = 1:params.Nms
% 
%             %Calculate current mean diameters
%             params.V_mu(iz, ix) = (params.nms(ix) .* (1 + params.X_mu(iz,ix)) .* params.R .* params.T_mu(iz,ix))./params.p_func(z);
%             params.d_mu(iz, ix) = ((6 .* params.V_mu(iz, ix))./pi).^(1/3);
% 
%             %Calculate current bubble velocities
%             params.uz_mu(iz, ix) = params.ubs.funcs{iz}(params.d_mu(iz, ix));
% 
%             %Calculate gas flow rates
%             V_dot(ix) = Ns_cell(ix) * V_cell * params.V_mu(iz, ix) * (params.uz_mu(iz,ix)./params.mesh.hs(iz));
%             V_gas = V_gas + Ns_cell(ix) * V_cell * params.V_mu(iz, ix);
% 
%         end
% 
%         V_dot_total = sum(V_dot); %m^3/s
% 
%         params.alpha_g(iz) = V_gas/V_cell;
% 
%         %Calculate gas superficial velocity in each cell
%         params.ug(iz) = V_dot_total/params.reactor.Ac;
% 
% 
%         %Calculate turbulent energy dissipation
%         params.turb.eps(iz) = params.g * params.ug(iz);
%         params.turb.k(iz) = 0.1;
%             %turb.eps_G(iz) = 1;
% 
%         %Check if value is reasonable
%         if params.alpha_g(iz) > params.alpha_g_max
%             params.alpha_g(iz) = params.alpha_g_max - 1E-6;
%             x = 1;
%         end
% 
%     end
% 
% end
% 
% 
% 
% function [c_src, c_snk] = Coalescence(y, params )
% 
%     %Check for inconsistent values
%     if any(y < 0)
%         warning('Negative numeric densities.')
%         y(y < 0) = 0;
%     end
% 
%     %Allocate vectors for source and sink temrs
%     c_src = zeros(size(y)); %Consider allocating empty vector in the params to avoid having to make a new one each iteration
%     c_snk = c_src;
% 
% 
%     %Iterate through all spatial cells
%     for iz = 1:params.Nz
% 
%         %Pull numeric densities in this spatial celll
%         cellinds = ((iz-1).*params.Nms + 1):1:(iz.*params.Nms);
%         Ns = y();
%         Ns_cell = Ns(cellinds);
% 
%         %Pull local fluid properties
%         uL = 0;
% 
% 
% 
% 
%         %Calculate average gas velocity from mass dependent velocities
%         % V_cell = params.Vcells(iz);
%         % V_gas = 0;
%         % for ix = 1:params.Nms
%         %     V_gas = V_gas + Ns_cell(ix) * V_cell * params.Vms(ix);
%         % 
%         % end
%         % params.alpha_g = V_gas/V_cell;
% 
%         %Calculate local fluid properties - only once per cell
%         turb = params.turb;
% 
% 
%         %Storage vector for sink term calculations
%         c_ij = zeros(params.Nms);
% 
% 
%         %Calculate all coalescence rates - needed for both sources/sinks
%         for j = 1:params.Nms
%             for k = 1:params.Nms
% 
%                 %Pull indices
%                 %i = xind;
%                 %j = coalesce_partners(ic,1);
%                 %k = coalesce_partners(ic,2);
% 
%                 %Caclulate Coalescence Rate
%                 c_ij(j,k) = Coalescence_Rate(Ns_cell,j,k, iz, params, turb);
%             end
%         end
% 
%         if any(c_ij(:) < 0)
%             x = 1;
%         end
% 
% 
%         %Iterate through internal variables within the cell
%         for ix = 1:params.Nms
% 
%             %Calculate indices
%             zind = iz;
%             xind = ix; %x refers to bubble size (mass currently)
%             yind = (zind-1) * params.Nms + xind;
% 
%             %Local density of this mass
%             N = y(yind); 
% 
%             Nb = N .* params.Vsz(zind);
% 
%             %Pull coalescence partners and iterate through them
%             coalesce_partners = params.coalesce_partners{xind};
%             coalesce_bias = params.coalesce_bias{xind};
%             coalesce_etas = params.coalesce_etas{xind};
%             coalesce_rats = params.coalesce_rats{xind};
% 
%             %Handle largest vin
%             if ix == params.Nms
%                x = 1; 
%             end
% 
%             %Iterate through 
%             cadd = 0;
%             csub = 0;
%             [N_partners, ~] = size(coalesce_partners);
%             for ic = 1:N_partners
% 
%                 %Pull indices
%                 i = xind;
%                 j = coalesce_partners(ic,1);
%                 k = coalesce_partners(ic,2);
%                 bias = coalesce_bias(ic);
%                 eta = coalesce_etas(ic);
% 
%                 %Pull original numeric densities
%                 jind = (zind-1) * params.Nms + j;
%                 kind = (zind-1) * params.Nms + k;
%                 Nj = y(jind);
%                 Nk = y(kind);
% 
%                 %Caclulate Coalescence Rate - pull from matrix
%                 c_jk = c_ij(j,k);   %Coalescence_Rate(y,i,j,k, xind, zind, params, turb);
% 
%                 %Determine rate for each of the three cells
%                 if bias < 0 %Shares with smaller (left) brackt
%                     c_jk_smaller = (1-eta) * c_jk * Nj * Nk;
%                     c_jk_middle = eta * c_jk * Nj * Nk;
%                     c_jk_larger = 0;
% 
%                 else %Shares with larger (right) bracket
%                     c_jk_smaller = 0;
%                     c_jk_middle = eta * c_jk * Nj * Nk;
%                     c_jk_larger = (1-eta) * c_jk * Nj * Nk;
%                 end
% 
%                 %Store results
%                 c_src(yind) = c_src(yind) + c_jk_middle;
%                 if ix > 1
%                     c_src(yind-1) = c_src(yind-1) + c_jk_smaller;
%                 end
%                 if ix < params.Nms
%                     c_src(yind+1) = c_src(yind+1) + c_jk_larger;
%                 end
% 
%                 %Subtract bubbles from their original bins
%                 c_snk(jind) = c_snk(jind) + c_jk * Nj * Nk;
%                 c_snk(kind) = c_snk(kind) + c_jk * Nj * Nk;
% 
%                 %
% 
% 
% 
%                 %Calculate total change for this bracket
%                % c_src(i) = c_src(i) + c_jk;
% 
%                 %Log result for sink calculations
%                % c_ij(j,k) = c_jk;
% 
%             end
% 
% 
% 
%         end
% 
%         %c_snk(i) = c_snk(i) - 0; %sum(Ns .* );
% 
% 
% 
%     end
% 
%     %Catch negative values
%     if any(c_src < 0)
%         x = 1;
%     end
% 
% end
% 
% 
% function c_jk = Coalescence_Rate(y, j,k, zind, params, turb)
% 
%     %Define index to sample from for liquid properties
%     i = (zind - 1) * params.Nms + 1;
% 
%     %Pull mass values at pivots
%         %mi = params.mms(i);
%     mj = params.mms(j);
%     mk = params.mms(k);
% 
%     %Determine which bracket 
% 
%     %Pull other values at pivots
%     jind = j; %jind = (zind - 1) .* params.Nms + j;
%     kind = k; %(zind - 1) .* params.Nms + k;
%     Nj = y(jind);
%     Nk = y(kind);
%     %Dd_j = ((6 .* params.mbs(j+1) .* params.R .* params.Ts(i))./(pi .* params.ps(i) .* params.Ms(i))).^(1/3)...
%       %  - ((6 .* params.mbs(j) .* params.R .* params.Ts(i))./(pi .* params.ps(i) .* params.Ms(i))).^(1/3);
%     %Dd_k = ((6 .* params.mbs(k+1) .* params.R .* params.Ts(i))./(pi .* params.ps(i) .* params.Ms(i))).^(1/3)... 
%       %  - ((6 .* params.mbs(k) .* params.R .* params.Ts(i))./(pi .* params.ps(i) .* params.Ms(i))).^(1/3);
% 
%                 % %Pull mass values of neighboring brackets
%                 % mi_lo = params.mms(i-1);
%                 % if i == params.Nms
%                 %     mi_hi = params.mbs(i+1);
%                 % else
%                 %     mi_hi = params.mms(i+1);
%                 % end
% 
%     %Calculate coalesced mass
%     mjk = mj + mk;
% 
%     %Calculate diameters of bubbles 
%         %di = ((6 * mi * params.R .* params.Ts(i))./(params.ps(i) .* params.Ms(i))).^(1/3);
%     dj = ((6 * mj * params.R .* params.Ts(i))./(pi .* params.ps(i) .* params.Ms(i))).^(1/3);
%     dk = ((6 * mk * params.R .* params.Ts(i))./(pi .* params.ps(i) .* params.Ms(i))).^(1/3);
% 
%     %Coalescence due to turbulent eddies - Adapted from Wang et al. 2005 -
%     %   Theoreetical prediction of flow regime ...
%     u_bar_j = sqrt(2) .* (turb.eps(zind) .* dj).^(1/3); %Mean turbulent velocity 
%     u_bar_k = sqrt(2) .* (turb.eps(zind) .* dk).^(1/3);
%     nj = Nj/params.mds(j);
%     nk = Nk/params.mds(k);  
%     hb_jk = 6.3 .* (nj + nk).^(-1/3); 
%     lbt_jk = sqrt( (0.89 .* dj).^2 + (0.89 .* dk).^2);
%     gamma_jk = exp(-(hb_jk./lbt_jk).^6);
%     omega_c = (pi./4) .* ((params.alpha_g_max)./(params.alpha_g_max - params.alpha_g(zind))) .* gamma_jk...
%         .* sqrt(2) .* turb.eps(zind).^(1/3) .* (dj + dk).^2 .* sqrt(dj.^(2/3) + dk.^(2/3));  %sqrt(u_bar_j.^2 + u_bar_k.^2);
%     psi = 1;
%     u_bar_jk = sqrt(u_bar_j.^2 + u_bar_k.^2);
%     We_ij = (params.rhos_l(i) .* (dj) .* u_bar_jk.^2)./params.sigmas(zind);
%     rho_g = 0.5;
%     gamma_VM = 0.5; %coefficient of virtual mass
%     Pc_jk = exp(-psi .* sqrt(0.75 .* (1+ (dj/dk).^2) .* (1 + (dj./dk).^3))./((rho_g/params.rhos_l(i) + gamma_VM) .* (1 + (dj./dk)).^3) .* sqrt(We_ij));
%     cc_jk = omega_c * Pc_jk;
% 
% 
%     %Coalescence due to wake entrainment
%     K_w1 = 15.4;
%     K_w2 = 0.46;
%     dc_w = 4 .* sqrt(params.sigmas(zind)./(params.g .* params.rhos_l(i))); %Critica diameter for wake entrainment
%     if dj  < dk
% 
%         %Calculate slip velocity - seems to be for a cap
%         %bubble
%         u_bar_slip = 0.71 .* sqrt(params.g.*dk);
% 
%         %Determine if bubble is above critical size for wake
%         %entrainment
%         if dk >= dc_w/2
%             theta_c = ((dk - dc_w/2).^6)./((dk - dc_w/2).^6 + (dc_w/2).^6);
%         else
%             theta_c = 0;
%         end
% 
%         %Calculate coalescence rate
%        % omega_w = K .* dk.^2 .* u_bar_slip;
%         cw_jk = K_w1 .* theta_c .* dk.^2 .* u_bar_slip .* exp(-K_w2 .* (((params.rhos_l(i).^(3) .* turb.eps(zind).^2)./(params.sigmas(zind).^3)) .* ((dj .* dk)./(dj + dk)).^5).^(1/6));
% 
%     else % j is the larger or equivalent bubble
% 
%         %Calculate slip velocity
%         u_bar_slip = 0.71 .* sqrt(params.g.*dj);
% 
%         %Determine if bubble is above critical size for wake
%         %entrainment
%         if dj >= dc_w/2
%             theta_c = ((dj - dc_w/2).^6)./((dj - dc_w/2).^6 + (dc_w/2).^6);
%         else
%             theta_c = 0;
%         end
% 
%         %Calculate coalescence rate
%        % omega_w = K_ .* dj.^2 .* u_bar_slip;
%         cw_jk = K_w1 .* theta_c .* dj.^2 .* u_bar_slip .* exp(-K_w2 .* (((params.rhos_l(i).^(3) .* turb.eps(zind).^2)./(params.sigmas(zind).^3)) .* ((dk .* dj)./(dk + dj)).^5).^(1/6));
%     end
% 
%     %Calculate probability for bubbles to coalesce 
% 
% 
%     %Calculate coalescence rate total
%     cjk = 0;
%     if params.coalesce.eddy
%         cjk = cjk + cc_jk;
%     end
%     if params.coalesce.wake
%         cjk = cjk + cw_jk;
%     end
%     %cjk = cc_jk + cw_jk; % Coalescence rate #/m^3-s - Used direclty for coalescence sink term
% 
%     %Calculate distribution to discrete brackets
%     if j == k
%         dcoeff = (1-0.5);
%     else
%         dcoeff = 1;
%     end
% 
%     %Pull distribution coeff
%     eta_ijk = 1;
% 
%     %Calculate distribution coeff - MOVE OUTSIDE INTEGRATION
%     %LOOP - PUT IN PBM
%     % if mjk < mi
%     %     eta_ijk = (mjk - mi_lo)./(mi - mi_lo);
%     % elseif mjk >= mi
%     %     eta_ijk = (mi_hi - mjk)./(mi_hi - mi);
%     % end
% 
%     %Calculate total change for this bracket
%     %c_src(i) = c_src(i) + dcoeff * eta_ijk * cjk * y(jind) * y(kind);
% 
% 
%     c_jk = cjk;
% 
% 
%     %c_jk = dcoeff * eta_ijk * cjk * y(jind) * y(kind); %USE FOR
%     %DISTRIBUTIN
% 
%       %  (1 - 0.5 * )
%     %eta_ijk = 1;
% 
%     if c_jk < 0
%         x = 1;
%     end
% 
% 
% 
% end
% 
% 
% function cc_jk = Coalescence_Turb(eps_turb, dj, dk)
% 
%     u_bar_j = sqrt(2) .* (eps_turb .* dj).^(1/3);
%     u_bar_k = sqrt(2) .* (eps_turb .* dk).^(1/3);
%     hb_jk = 6.3 .* (Nj/Dd_j + Nk/Dd_k).^(-1/3);
%     lbt_jk = sqrt( (0.89 .* dj).^2 + (0.89 .* dk).^2);
%     gamma_jk = exp(-(hb_jk./lbt_jk).^6);
%     omega_c = (pi./4) .* ((alpha_g_max)./(alpha_g_max - alpha_g)) .* gamma_jk...
%         .* sqrt(2) .* eps_turb.^(1/3) .* (dj + dk).^2 .* sqrt(dj.^(2/3) + dk.^(2/3));  %sqrt(u_bar_j.^2 + u_bar_k.^2);
%     psi = 1;
%     u_bar_jk = sqrt(u_bar_j.^2 + u_bar_k.^2);
%     We_ij = (params.rhos_l(i) .* ((dj + dk)./2) .* u_bar_jk.^2)./params.sigmas(i);
%     rho_g = 0.5;
%     Pc_jk = exp(-psi .* sqrt(0.75 .* (1+ (dj/dk).^2) .* (1 + (dk./dk).^3))./((rho_g/params.rhos_l(i) + 0.5) .* (1 + (dj./dk)).^3) .* sqrt(We_ij));
%     cc_jk = omega_c * Pc_jk;
% 
% 
% end
% 
% 
% function Coalescence_Wake
% 
% 
% 
% end
% 
% 
% function [b_src, b_snk] = Breakage(y, params)
% 
%     %Debug message
%     if params.debug
%         t_start = cputime;
%         fprintf('--Breakage Start\n');
%     end
% 
%     %Settings
%     params.delta = 1E-3; %Value of 0.01 is used in Wnag et al (2003)
% 
% 
%     %debug dettings
%     if params.debug
%         %params.turb.eps = 1 .* ones(size(params.turb.eps));
%     end
% 
%     %Pull values and allocate output vectors
%     b_src = zeros(size(y)); b_snk = b_src;
% 
% 
%     %Calculate kolmogorov length scale
%     lambda_komogorov = ((params.nus.^3)./params.turb.eps).^0.25; %m
%     lambda_min = 31.4 * lambda_komogorov; %Minimum eddy diameter to break a bubble - integrate up to bubble diameter 
% 
% 
%     %Iterate through all spatial cells
%     for iz = 1:params.Nz
% 
%         %Pull numeric densities in this spatial celll
%         z = params.zms(iz);
%         cellinds = ((iz-1).*params.Nms + 1):1:(iz.*params.Nms);
%         Ns = y();
%         Ns_cell = Ns(cellinds);
% 
%         %Pull local fluid properties
%         uL = 0;
% 
%         %Calculate critical breakup diameter
%         d_crit = CalcCritDiameter(iz, params); %m - critical diameter for instability breakage
% 
%         %Create storage matrix
%         bs = zeros(params.Nms, params.Nms); %Row = bubble breaking, column = bubble receiving
% 
% 
%         %VARY RESOLUTION OF fvs norm based on the number of size groups
%         %below it to reduce computational cost
% 
%         %Iterate through size groups
%         betas = zeros(params.Nms-1, length(params.fvs_norm_all));
%         for im = 2:params.Nms %params.Nms %Neglect first size group, since no where for it to break
% 
%             %Pull index
%             mind = cellinds(im);
% 
%             %Pull bubble size
%             mi = params.mms(im); %Current represenative mass
%             ni = params.nms(im);
%             V = (ni .* (1 + params.X_mu(iz, im)) .* params.R .* params.T_mu(iz,im))./params.p_func(z);
%             d = (6.*V/pi).^(1/3);
%             u = params.ubs.funcs{iz}(d); %params.uzs(cellinds(im));
% 
%             %Eddy Breakage
%             if params.break.eddy
%                 if params.break.interp %Solve using interpolation
%                     if any(~isreal(params.turb.eps(iz)))
%                         x = 1;
%                     end
% 
% 
%                     beta_ratio = params.break.beta_ratio{iz}(d, params.turb.eps(iz));
%                         %b_eddy = params.break.b_eddy{iz}(d, params.turb.eps(iz));
%                     beta_eddy = params.break.beta{iz}(d, params.turb.eps(iz), params.fvs_norm_all);
%                     beta_eddy = beta_eddy(:);
%                     if any(isnan(beta_eddy))
%                         beta_eddy = zeros(size(beta_eddy));
%                         beta_ratio = 0;
%                     end
%                     b_eddy = BreakageEddySimple(iz, im, d, Ns_cell, beta_ratio, beta_eddy, params);
%                 else %Solve in real time
%                     [b_eddy, beta_eddy] = BreakageEddyAlt(iz, im, d, Ns_cell, lambda_min, params.N_lambdas, params); %b_eddy = breakage rate of the size; beta_eddy = distribution of bubble sizes
%                 end
%             end
% 
%             if strcmp(params.break.type, 'uniform')
%                 beta_eddy_orig = beta_eddy;
%                 beta_eddy = ones(size(beta_eddy));
%             end
% 
%             %Inertia Breakage - Bubble too large
%             if params.break.surf
%                 if d > d_crit
%                     [b_surf, beta_surf] = BreakageSurf(d, d_crit, params);
%                 else
%                     b_surf = 0;
%                     beta_surf = 0;
%                 end
%             else
%                 b_surf = 0;
%                 beta_surf = 0;
%             end
% 
%             %Log BSD
%             betas(im-1, :) = beta_eddy;
% 
%             %SOURCE TERM --------------------------------------------------
%             b_surf = b_surf .* Ns_cell(im);
%             b_total = b_eddy + b_surf;
% 
% 
%             %Only consider if breakage rate is above a threshold
%             if b_total > 1E-6
% 
% 
% 
%                 %Calculate distribution of products
%                 ms_norm = mi * [params.fvs_norm, 1 - fliplr(params.fvs_norm(1:end-1))];
%                 ms_max = max(ms_norm);
%                     %ds_norm = (d.^3 .* params.fvs_norm).^(1/3); %Can calculate and list as a matrix with N_fvs columns, and Nms rows 
%                 zetas = zeros(1,im);
% 
%                 %Iterate through smaller brackets
%                 for is = 1:im %Only consider sizes smaller than 
% 
% 
%                     %SEARCH HERE FOR MISTAKE
% 
%                     %Pull cellind
%                     sind = cellinds(is); %source index
% 
%                     %Create interpolation function
%                     beta_func = @(m) (1/ms_max) .* interp1(ms_norm, beta_eddy, m); %Normalized to mass 
% 
%                     %Specify representative masses of adjacent cells (with
%                     %   which the new population is shared)
%                     if is > 1
%                         m_low = params.mms(is-1);
%                     else
%                         m_low = 0;
%                     end
%                     if is < params.Nms
%                         m_hig = params.mms(is+1);
%                     end
%                     m_mid = params.mms(is);
% 
%                     %Calculate integral
%                     %if im == 1 && is == 1
% 
%                     if is == 1 %Only preserve mass, not numbers
%                         int_func_up = @(m) (m_hig - m)./(m_hig - m_mid) .* beta_func(m);
%                         int_func_down = @(m) (m - 0)./(m_mid) .* beta_func(m); %Will conserve mass, but not number of bubbles since there is no lower bracket
%                         zetas(is) = integral(int_func_up, m_mid, m_hig) + integral(int_func_down, 0, m_mid);
%                     elseif is == params.Nms
%                         x = 1;
%                     elseif is == im
%                         int_func_down =  @(m) (m - m_low)./(m_mid - m_low) .* beta_func(m);
%                         zetas(is) = integral(int_func_down, m_low, m_mid); %Don't consider up direction since bubbles cannot be larger;
%                     else
%                         int_func_up = @(m) (m_hig - m)./(m_hig - m_mid) .* beta_func(m);
%                         int_func_down = @(m) (m - m_low)./(m_mid - m_low) .* beta_func(m);
%                         zetas(is) = integral(int_func_up, params.mms(is), params.mms(is+1)) + integral(int_func_down, params.mms(is-1), params.mms(is));
%                     end
% 
%                     %Add to source term
%                     b_src(sind) = b_src(sind) + zetas(is) .* b_eddy + 1/(im-1) * b_surf;
% 
% 
% 
% 
%                 end
% 
% 
%             end
% 
%             %SINK TERM ----------------------------------------------------
%             b_snk(mind) = b_total;
% 
% 
% 
%         end
% 
%         %Debug plots
%         if params.break.debug & false
% 
%             %Create initial figure
%             figure('units', 'normalized', 'OuterPosition',  [0, 0, 1, 1]);
%             linespecs = {'-', '--', ':', '-.'};
%             it = 2;
%             for i = 2:5; %params.Nms
% 
%                 %Plot BSD
%                 %subplot(3,5,it);
%                 plot(params.fvs_norm_all, betas(i-1,:), 'k-', 'LineWidth', 1.5, 'LineStyle', linespecs{i-1}); hold on;
%                 grid on; grid minor; axis square;
%                 xlabel('fv - breakup vol. fraction'); ylabel('Normalized BSD');
%                 title(sprintf('BSD (im = %d)', i));
%                 it = it + 1;
% 
%                 %Create new figure and reset iteration counter
%                 if it > 15
%                     figure('units', 'normalized', 'OuterPosition',  [0, 0, 1, 1]);
%                     it = 1;
%                 end
%             end
% 
%             %Plots
%             title('Normal Bubble Size Distribution');
%             legend('d = 1.5 mm', 'd = 2 mm', 'd = 3 mm', 'd = 6 mm');
%             set(gca, 'FontSize', 18);
% 
% 
%         end
% 
%     end
% 
%     %Debug message
%     if params.debug
%         t_end = cputime;
%         t_req = t_end - t_start;
%         fprintf('--Breakage End (t = %0.4f s)\n', t_req);
%     end
% 
% 
% 
% 
% end
% 
% 
% 
% function [b_surf, fv] = BreakageSurf(d, d_crit, params)
% 
%     %Additional Parameters - pronably move outside loop
% 
% 
%     %
%     b_surf = params.break.b_star .* ((d - d_crit).^params.break.m_star)./((d - d_crit).^params.break.m_star + d_crit.^params.break.m_star);
%     fv = 0.5; %0.5 = equal breakup
% 
% end
% 
% 
% function [b_eddy, beta] = BreakageEddy(iz, im, d, Ns_cell, lambda_min, N_lambdas, params)
% 
%     %Debug cases
%     if params.break.debug
%         %d = 0.003;
%         t_start.total = cputime;
%         fprintf('-Eddy Breakage: iz = %i; im = %i; d = %0.5f m; m_i = %0.3e;\n', iz, im, d, params.mms(im));
%         x = 1;
%     end
% 
%     %Specify standard fv distribution to interpolate onto
%     % fvs_norm = [0:params.dfv_cap:0.05, (0.05+params.dfv_surf):params.dfv_surf:0.5];
% 
%     %Create storage vectors
%     Pbs_norm = params.Pbs_norm; %zeros(N_lambdas, length(params.fvs_norm));
% 
%     %Iterate through lambdas
%     it_total = 0;
%     if lambda_min(iz) < d %otherwise the bubble cannot be broken by eddies
% 
%         %
%         lambdas = linspace(lambda_min(iz), d, N_lambdas);
% 
%         %Calculate cutoffs for each lambda
%         u_bar_lambdas = sqrt(2) .* (params.turb.eps(iz) .* lambdas).^(1/3);
%         e_bar_lambdas = (pi/6) .* lambdas.^3 .* params.rhos(iz) .* (u_bar_lambdas.^2)./2; %Might be able to move this calculation outside of ODE func - use interp
%         e_cutoffs = 10 .* e_bar_lambdas;
% 
%         for il = 1:N_lambdas
% 
%             %Debug timing
%             if params.debug
%                 t_start.setup = cputime;
%             end
% 
%             %Pull current eddy diamter
%             lambda = lambdas(il);
% 
%             lambda_rat =  lambda/d;
%             %lambda = lambda_rat * d;
% 
%             %Evaluate transition break fraction 
%             fvc = params.fvc_func(lambda_rat);
% 
%             %Specify values to integrate over
%             fvs_surf = 0.5:(-params.dfv_surf):fvc; fvs_surf(end+1) = fvc;
%             fvsm_surf = (fvs_surf(1:end-1) + fvs_surf(2:end))./2;
%             if fvc/params.dfv_cap > params.N_fv_min
%                 fvs_cap = 1E-8:params.dfv_cap:fvc;
%             else
%                 fvs_cap = linspace(1E-8,fvc,params.N_fv_min);
%             end
%             fvsm_cap = (fvs_cap(1:end-1) + fvs_cap(2:end))./2;
% 
%             %Calculate average velocity and energy of eddy
%             %u_bar_lambda = sqrt(2) .* (params.turb.eps(iz) .* lambda).^(1/3);
%             %e_bar_lambda = (pi/6) .* lambda.^3 .* params.rhos(iz) .* (u_bar_lambda.^2)./2;
%             %e_cutoff = 10 * e_bar_lambda;
% 
%             %Calculate min and max breakup fraction
%             cf_max = min([(2^(1/3) - 1), e_cutoffs(il)./(pi .* d.^2 .* params.sigmas(iz))]);
%             fv_func = @(fv) fv.^(2/3) + (1-fv).^(2/3) - 1 - cf_max;
%             fv_max = fzero(fv_func, [0, 0.5]);
%             fv_max = round(fv_max,3);
% 
%             if fv_max > 0.49
%                 x = 1;
%             end
% 
%             % fv_min = ((pi .* lambda.^3.*params.sigmas(iz))./(6 * e_cutoff * d)).^3;
%             % cf_max = min([(2.^(1/3) - 1), e_cutoff./(pi .* d.^2 .* params.sigmas(iz))]);
%             % fv_func = @(fv) fv.^(2/3) + (1-fv).^(2/3) - 1 - cf_max;
%             % fv_max = fzero(fv_func, 0.5);
%             % cf_half = 2^(1/3) - 1;
%             % e_half = cf_half * pi * d^2 * params.sigmas(iz);
% 
%             %Calculate energy required to break the bubble into
%             %each possible fraction
%             cfs = (fvs_surf.^(2/3) + (1-fvs_surf).^(2/3) - 1);
%             e_is = cfs .* pi .* d.^2 .* params.sigmas(iz); %Energy required to break - should be calculated outside of loop
%             e_ihs = (e_is(1:end-1) + e_is(2:end))./2;
% 
%             %Calculate initial value at 0.5
%             if e_is(1) < e_cutoffs(il)
%                 fv = 0.5;
%                 e_min = (fv.^(2/3) + (1-fv).^(2/3) - 1) .* pi .* d.^2 .* params.sigmas(iz);
%                 e_max = e_cutoffs(il);
%                 fv_min = ((pi .* lambda.^3 .* params.sigmas(iz))./(6 .* e_ihs(1) .* d)).^3;
%                 Pbe = 1./(fv - fv_min);
%                 Pe = (1./e_bar_lambdas(il)) .* exp(-e_ihs(1)/e_bar_lambdas(il)); %
%                 Pb_surf_init = Pbe .* Pe .* (e_is(1) - e_is(2));
% 
% 
%                 P_func = @(e) 1./(fv - fv_min) .* (1./e_bar_lambdas(il)) .* exp(-e/e_bar_lambdas(il));
% 
%                 Pb_surf_init = integral(P_func, e_is(1), e_max);
% 
%             else
%                 Pb_surf_init = 0;
%             end
% 
% 
%              %To calculate initial Pb at fv = 0.5 for the surface controlled region, need to consider cases fv_max = 0.5. If fv_max < 0.5, then it starts at zero as it currently does, otherwise the baseline is higher.
% 
%             %Debug timing
%             if params.debug
%                 t_end = cputime;
%                 t_req.setup = t_end - t_start.setup;
%                 fprintf('--Eddy Breakage | Setup time   = %0.6f s\n', t_req.setup);
%                 t_start.surf = cputime;
%             end
% 
%             %Iterate through possible breakup fractions
%             Pb_surf = zeros(1, length(fvsm_surf));
%             can_break_surf = e_ihs <  e_cutoffs(il); %Logical stating if any eddies have the energy to break the bubble into the defined fractions
% 
%             %Surface energy region ----------------------------------------
%             if any(can_break_surf) %Check if no eddies can break this bubble
% 
% 
% 
%                 %Iterate only through cases which can cause
%                 %breakage
%                 break_inds = find(can_break_surf);
%                 for iv = break_inds
% 
%                     %Calculate probability of eddy having the required kinetic energy
%                     Pe = (1./e_bar_lambdas(il)) .* exp(-e_ihs(iv)/e_bar_lambdas(il)); %probability of an eddy with this amount of energy
% 
%                         %Vectorize and do outside loop
% 
%                     fv = fvs_surf(iv);
%                     %fvh = (fvs_surf(iv) + fvs_surf(iv+1))./2;
%                     %cf = (fv.^(2/3) + (1-fv).^(2/3) - 1);
%                     %efv = cf * pi * d^2 * params.sigmas(iz);
% 
% 
%                     %Calculate minimum energy and minimum breakup fraction
%                     e_min = (fv.^(2/3) + (1-fv).^(2/3) - 1) .* pi .* d.^2 .* params.sigmas(iz);
%                     fv_ih_min = ((pi .* lambda.^3 .* params.sigmas(iz))./(6 .* e_ihs(iv) .* d)).^3;
% 
%                     %Calculate breakage probability
%                     if (fvsm_surf(iv) - fv_ih_min) >= params.break.delta
%                         Pbe = 1/(fvs_surf(iv) - fv_ih_min);
%                     else
%                         Pbe = 0;
%                     end
% 
%                     %Calculate total probability
%                     if iv == 1
%                         Pb_surf(iv) = (e_is(iv) - e_is(iv+1)) .* Pbe .* Pe + Pb_surf_init;
%                     else
%                         %Pb_surf(iv) = (e_is(iv) - e_is(iv+1)) .* Pbe .* Pe + Pb_surf(iv-1);
% 
%                         %Alternative calculation
%                         P_func = @(e) 1./(fv - fv_ih_min) .* (1./e_bar_lambdas(il)) .* exp(-e/e_bar_lambdas(il));
%                         Pb_surf(iv) = integral(P_func, e_ihs(iv), e_cutoffs(il)) + Pb_surf(iv-1);
%                     end  
% 
% 
% 
%                     %Log number of iterations
%                     it_total = it_total + 1;
% 
%                 end
% 
%                 %Fill rest of 
%             end
% 
%             %Debug timing
%             if params.debug
%                 t_end = cputime;
%                 t_req.surf= t_end - t_start.surf;
%                 fprintf('--Eddy Breakage | Surf. time   = %0.6f s\n', t_req.surf);
%                 t_start.cap = cputime;
%             end
% 
% 
% 
% 
%             %Pressure controlled region --------------------------
%             e_is = (pi .* lambda.^3 .* params.sigmas(iz))./(6 .* fvs_cap.^(1/3) .* d);
%             e_ihs = (e_is(1:end-1) + e_is(2:end))./2;
% 
%             %Iterate through possible breakup fractions
%             Pb_cap = zeros(1, length(fvsm_cap));
%             can_break_cap = e_ihs < e_cutoffs(il); 
%             if any(can_break_cap)
%                 break_inds = find(can_break_cap);
%                 for iv = break_inds
% 
%                     %Calculate Pbe - Pb for the specific energy
%                     fv_ih_max_func = @(fv) e_ihs(iv) - (fv.^(2/3) + (1 - fv).^(2/3) - 1) .* pi .* d.^2 .* params.sigmas(iz);
%                     try
%                         fv_ih_max = fzero(fv_ih_max_func, [0, 0.5]);
%                         %fv_ih_max = params.fvih_max_funcs{iz}(d, e_ihs(iv)); %fzero(fv_ih_max_func, 0.05);
%                         if fv_ih_max > 0 && (fv_ih_max - fvsm_cap(iv)) >= params.delta
%                             Pbe = 1/(fv_ih_max - fvsm_cap(iv));
%                         else
%                             Pbe = 0;
%                         end
%                     catch
%                         Pbe = 0;
%                     end
% 
%                     %Calculate 
%                     Pe = (1./e_bar_lambdas(il)) .* exp(-e_ihs(iv)/e_bar_lambdas(il));
% 
%                      %Calculate total probability
%                     if iv == 1
%                         Pb_cap(iv) = (e_is(iv) - e_is(iv+1)) .* Pbe .* Pe;
%                     else
%                         Pb_cap(iv) = (e_is(iv) - e_is(iv+1)) .* Pbe .* Pe + Pb_cap(iv-1);
%                     end  
% 
% 
%                     if any(Pb_cap < 0)
%                         x = 1;
%                     end
% 
%                     %Log number of iterations
%                     it_total = it_total + 1;
% 
%                 end
% 
%             end
% 
%             %Debug timing
%             if params.debug
%                 t_end = cputime;
%                 t_req.cap= t_end - t_start.cap;
%                 fprintf('--Eddy Breakage | Cap. time    = %0.6f s\n', t_req.cap);
%                 t_start.interp = cputime;
%             end
% 
% 
% 
%             %Interpolate results
%             fvs = [fvsm_cap, fliplr(fvsm_surf)];
%             [fvs, unique_inds] = unique(fvs);
%             Pbs = [Pb_cap, fliplr(Pb_surf)];
%             Pbs = Pbs(unique_inds);
%             Pbs_norm(il,:) = interp1(fvs, Pbs, params.fvs_norm, 'linear', 'extrap');
%             Pbs_norm(il,1) = 0;
%             % 
%             %Debug analysis
%             if params.debug
% 
%                 %Plot
%                 figure();
%                 plot(params.fvs_norm, Pbs_norm(il,:), 'k-', 'LineWidth', 2);
%                 grid on; grid minor; axis square;
% 
%                 if lambda_rat > 0.69
%                     x =1;
%                 end
%                 close 
% 
%             end
% 
%             %Debug timing
%             if params.debug
%                 t_end = cputime;
%                 t_req.interp = t_end - t_start.interp;
%                 fprintf('--Eddy Breakage | Interp. time = %0.6f s\n', t_req.interp);
%             end
% 
% 
% 
%             %Calculate complete integral
%             %int_lambda(il) = sum(Pbs_norm(il,:)) .* ((lambdas(il) + d).^2)./(lambdas(il).^(11/3)); 
% 
% 
% 
%             %  %Debug updates
%             % if params.debug
%             %     fprintf('im = %d; d = %0.4f mm; iv = %d; Pb = %0.4e;\n', im, d, iv, Pb)
%             % end
%         end
% 
%     else
% 
%         %Probability is zero - lambda is too large, will just transport
%         %bubble 
%         error('Finish this.');
% 
%     end
% 
% 
%     %Evaluate outer integral
%     b_fvd = params.bfd_zero;
%     for iv = 2:length(params.fvs_norm)
%         f_lambda_fv = Pbs_norm(:,iv)' .* ((lambdas + d).^2)./(lambdas.^(11/3));
%         int_lambda_fv = trapz(lambdas, f_lambda_fv);
%         b_fvd(iv) = 0.923 .* (1 - params.alpha_g(iz)) .* Ns_cell(im) .* params.turb.eps(iz).^(1/3) .* int_lambda_fv;
% 
%         %Log number of iterations
%         it_total = it_total + 1;
%     end
% 
%     %Calculate overall rate
%     b_eddy = trapz(params.fvs_norm, b_fvd);
% 
%     %Determine normalized daughter bubble size distribution
%     beta = (b_fvd)./(b_eddy);
%     beta = [beta, fliplr(beta(1:end-1))];
% 
%     %
%     if any(isnan(beta))
%         beta = zeros(size(beta));
%     end
% 
% 
%     %Debug plot
%     if params.break.debug
%         t_end = cputime;
%         t_req = t_end - t_start.total;
%         figure();
%         plot(params.fvs_norm_all, beta);
%         close
%         fprintf('-- Eddy Breakage time = %0.4f; Total its = %d; \n', t_req, it_total)
%     end
% 
% end
% 
% 
% %Calculate the breakage rate from interpolated values of 
% function b_eddy = BreakageEddySimple(iz, im, d, Ns_cell, beta_ratio, beta_eddy, params)
% 
%     %Debug cases
%     if params.break.debug
%         %d = 0.003;
%         t_start.total = cputime;
%         fprintf('-Eddy Breakage: iz = %i; im = %i; d = %0.5f m; m_i = %0.3e;\n', iz, im, d, params.mms(im));
%         x = 1;
%     end
% 
%     %Track iterations
%     it_total = 1;
% 
%     %Calculate integral
%     b_fvd = params.break.bfd_zero;
%     int_lambda_fvs = beta_ratio .* beta_eddy(1:length(params.fvs_norm));
%     b_fvd = 0.923 .* (1 - params.alpha_g(iz)) .* Ns_cell(im) .* params.turb.eps(iz).^(1/3)  .* int_lambda_fvs;
% 
%     % for iv = 2:length(params.fvs_norm)
%     % 
%     %     %Calulcate 
%     %     int_lambda_fv = beta_ratio .* beta_eddy(iv);
%     %     b_fvd(iv) = 0.923 .* (1 - params.alpha_g(iz)) .* Ns_cell(im) .* params.turb.eps(iz).^(1/3) .* int_lambda_fv;
%     % 
%     %     %Update iteration counter
%     %     it_total = it_total + 1;
%     % 
%     % end
% 
%     %Calculate overall rate
%     b_eddy = trapz(params.fvs_norm, b_fvd);
% 
%     %Debug plot
%     if params.break.debug
%         t_end = cputime;
%         t_req = t_end - t_start.total;
%         figure();
%         plot(params.fvs_norm_all, beta_eddy);
%         close
%         fprintf('-- Eddy Breakage time = %0.4f; Total its = %d; \n', t_req, it_total)
%    end
% 
% 
% 
% end
% 
% 
