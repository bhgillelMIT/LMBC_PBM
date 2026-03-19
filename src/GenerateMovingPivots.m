
function gparams = GenerateMovingPivots(gparams, params)

    %Handle no inpits
    if nargin < 1
        PBM_v3()
    end

    %Physical constants
    R = 8.3145;
    g = 9.81;

    %Pull values
    T_mu = gparams.T_mu;
    T_std = gparams.T_std;
    N_bs = gparams.N_bs;

    %Core Distribution parameters
    gparams.core.N_stds_bracket = gparams.core.N_stds./gparams.core.N_pivots;
    gparams.core.T_min = T_mu - gparams.core.N_stds * gparams.T_std;
    gparams.core.T_max = T_mu + gparams.core.N_stds * gparams.T_std;
    gparams.core.Tms = linspace(gparams.core.T_min, gparams.core.T_max, 2*gparams.core.N_pivots+1);
    gparams.core.dTm_init = mean(diff(gparams.core.Tms));
    gparams.core.Tps = (gparams.core.Tms(1) - gparams.core.dTm_init/2):gparams.core.dTm_init:(gparams.core.Tms(end) + gparams.core.dTm_init/2);
    gparams.core.Fps = N_bs/(T_std * sqrt(2*pi)) .* exp(-0.5 * ((gparams.core.Tms - gparams.T_mu).^2)./(gparams.T_std^2));
    gparams.core.T_min = min(gparams.core.Tps);
    gparams.core.T_max = max(gparams.core.Tps);
    

    %Add engaged trailing points
    gparams.core.dT_trail = gparams.core.T_min - gparams.T_min;
    gparams.core.N_trail_engaged = floor(gparams.core.dT_trail/gparams.trail.dT);
    gparams.core.Tps_trail = (gparams.core.T_min - gparams.core.N_trail_engaged * gparams.trail.dT):gparams.trail.dT:gparams.core.T_min;
    gparams.core.Tms_trail = gparams.core.Tps_trail + gparams.trail.dT/2;
    gparams.core.Tps = [gparams.core.Tps_trail, gparams.core.Tps];
    gparams.core.Tms = [gparams.core.Tms_trail, gparams.core.Tms];
    gparams.core.Fps_trail = zeros(size(gparams.core.Tms_trail));
    gparams.core.Fps = [gparams.core.Fps_trail, gparams.core.Fps];
    gparams.core.engaged = true(size(gparams.core.Tps));

    %Other parameters
    gparams.As = @(T,p,M) 4 * pi * ((3 * gparams.m .* R .* T)./(4.*pi.*p.*M)).^(2/3); %m^2 - Surface area as a function of temperature 
    gparams.A_func = @(T,p,M) (gparams.h(T) .* gparams.As(T,p,M))./(gparams.m .* gparams.Cp(T));  %@(T) (params.h .* params.As)./(gparams.ms(T) .* params.Cp);
    gparams.As_max = @(T,p,M) 4 * pi * ((3 * gparams.m_max .* R .* T)./(4.*pi.*p.*M)).^(2/3);
    gparams.A_min_func = @(T,p,M) (gparams.h(T) .* gparams.As_max(T,p,M))./(gparams.m_max .* gparams.Cp(T)); %Lowest heating rate coefficient - corresponds to largest bubble 
    params.gparams = gparams;

    %Determine trail points
    trail = DetermineMaxDomain(params); %Estimate velocity of slowest heating bubble
    gparams.trail.Tis = trail.Tis;
    gparams.trail.zis = trail.zis;
    gparams.trial.tis = trail.tis;
    gparams.trail.uis = trail.uis;
    gparams.trail.N_pivots = length(gparams.trail.Tis);
    gparams.trail.engaged = false(size(gparams.trail.Tis));
    gparams.trail.Fps = zeros(size(gparams.trail.engaged));
    gparams.trail.Tps = gparams.trail.Tis;
    gparams.trail.slowest = trail.slowest;





    % gparams.trail.dT_core = min(gparams.core.Tps) - gparams.T_min;
    % gparams.dT_core_max = gparams.trail.dT_core + gparams.dT_max;
    % gparams.trail.N_pivots = max([ceil(gparams.dT_core_max/gparams.trail.dT), 1]);
    % 
    % gparams.trail.N_engaged = floor(gparams.trail.dT_core/gparams.trail.dT);
    % 
    % 
    % gparams.trail.N_unengaged = gparams.trail.N_pivots - gparams.trail.N_engaged;
    % gparams.trail.Ts_min_engaged = gparams.core.T_min - gparams.trail.dT * gparams.trail.N_engaged;
    % gparams.trail.Ts_engaged = linspace(gparams.trail.Ts_min_engaged, gparams.core.T_min - gparams.trail.dT, gparams.trail.N_engaged);
    % gparams.trail.Ts_unengaged = gparams.T_min:gparams.trail.dT_init:(gparams.T_min + (gparams.trail.N_unengaged - 1)*gparams.trail.dT_init); %Unengaged points will travel at speed of slowest (i.e. largest) bubble
    % if ~isempty(gparams.trail.Ts_engaged)
    %     if any(gparams.trail.Ts_unengaged > min(gparams.trail.Ts_engaged))
    %         gparams.trail.dT_init = (min(gparams.trail.Ts_engaged) - gparams.T_min)/(gparams.trail.N_pivots - gparams.trail.N_engaged);
    %         gparams.trail.Ts_unengaged = gparams.T_min:gparams.trail.dT_init:(gparams.T_min + (gparams.trail.N_unengaged-1) * gparams.trail.dT_init);
    %     end
    % end
    % gparams.trail.Tps = [gparams.trail.Ts_unengaged, gparams.trail.Ts_engaged]; 
    % gparams.trail.Tms = (gparams.trail.Tps + [gparams.trail.Tps(2:end), gparams.core.T_min])./2;
    % gparams.trail.engaged = [false(1,gparams.trail.N_unengaged), true(1,gparams.trail.N_engaged)];
    % gparams.trail.Fps = N_bs/(T_std * sqrt(2*pi)) .* exp(-0.5 * ((gparams.trail.Tms - T_mu).^2)./(T_std^2));
        % gparams.Ts_trail = linspace(T_min, gparams.T_min, gparams.N_pivots_trail+1);
        % gparams.Ts_trail = gparams.Ts_trail(1:end-1);
        % gparams.Fs_trail = N_bs/(T_std * sqrt(2*pi)) .* exp(-0.5 * ((gparams.Ts_trail - T_mu).^2)./(T_std^2));

    %Determine lead points 
    %gparams.lead.dT = 200;
    gparams.lead.N_pivots = ceil((gparams.T_liq - max(gparams.core.Tps))/gparams.lead.dT);
    gparams.lead.Tps = gparams.T_liq:-gparams.lead.dT:gparams.core.T_max;
    gparams.lead.Tps = fliplr(gparams.lead.Tps(2:end));
    gparams.lead.Tms = ([gparams.core.T_max, gparams.lead.Tps(1:end-1)] + gparams.lead.Tps(1:end))./2;
    gparams.lead.Fps = N_bs/(T_std * sqrt(2*pi)) .* exp(-0.5 * ((gparams.lead.Tms - T_mu).^2)./(T_std^2));
                        %Last point (local liquid temperature) is implicit, doesn't need to be
                        %similated
    gparams.lead.engaged = true(size(gparams.lead.Fps));

    %Define all pivots
    gparams.Fps = [gparams.trail.Fps, gparams.core.Fps, gparams.lead.Fps];    
    gparams.Tps = [gparams.trail.Tps, gparams.core.Tps, gparams.lead.Tps];
    gparams.Tms = gparams.Tps(1:end-1) + diff(gparams.Tps)./2;
    gparams.engaged = [gparams.trail.engaged, gparams.core.engaged, gparams.lead.engaged];

    %Count
    gparams.N_brackets = length(gparams.Fps);
    gparams.N_Ts = length(gparams.Tps);
    gparams.N_Fs = length(gparams.Fps);

    %Calculate numeric densities
    Tms = [gparams.core.Tps, gparams.lead.Tps, gparams.T_liq];
    dist_func = @(T) 1./(T_std .* sqrt(2.*pi)) .* exp(-0.5 .* ((T - gparams.T_mu).^2)./(gparams.T_std^2));
    Fracs = zeros(1, length(Tms)-1);
    for it = 1:length(Fracs)
        Fracs(it) = integral(dist_func, Tms(it), Tms(it+1));
    end
    Fracs_trail = zeros(1,(length(gparams.trail.Tps)));
    gparams.Fracs = [Fracs_trail, Fracs];



end


%A function to identify the max separation between the 
function trail = DetermineMaxDomain(params)

    %Constant 


    %Physical constants
    R = 8.3145;

    T_i = params.T_min; %params.gparams.T_i;

    %Setup ode system
    odeparams.H_func = @(n_CH4, n_H2, n_C, T) n_CH4 * integral(methane.Cp, T0, T) + ...
        n_H2 * integral(hydrogen.Cp, T0, T) + n_C * integral(carbon.Cp, T0, T);
    odeparams.p_z = params.p_z;
    odeparams.h_reactor = params.reactor.H;
    odeparams.T_Lz = params.T_Lz;
    odeparams.rho_L_T = params.rho_L_T;
    odeparams.g = params.g;
    odeparams.R = params.R;
    odeparams.methane = params.methane;
    odeparams.hydrogen = params.hydrogen;
    odeparams.carbon = params.carbon;
    odeparams.liquud = params.liquid;
    odeparams.argon = params.argon;
    odeparams.k0 = 6.6E13;
    odeparams.Ea = -370000;
    odeparams.h_conv = 5;
    odeparams.Xis = 0:0.125:0.875;
    odeparams.chars = params.chars;

    
    %Calculate slowest characteristic for largest bubble
    m = params.mms(end);
    n_i = params.nms(end);
    n_Ar_i = n_i * params.X_Ar;
    n_CH4_i = n_i * (1-params.X_Ar);
    odeparams.h_conv_const = params.h_conv_const;
    odeparams.im = params.Nms;
    odeparams.m = m;
    odeparams.liquid = params.liquid;
    odeparams.Xi = 0; %
    odeparams.n_i = n_i;
    odeparams.n_Ar_i = n_Ar_i;
    odeparams.n_CH4_i = n_CH4_i;
    odeparams.m = params.mms(end);
    [T1,Y1] = TemperatureCharacteristic(0, params.t_f, T_i, [], params, odeparams, false);
    
    %Isolate values for largest bubble
    slowest.t = T1;
    slowest.u = Y1(:,1);
    slowest.z = Y1(:,2);
    slowest.T = Y1(:,3);
    slowest.X = Y1(:,4);
    slowest.T_end = max(slowest.T);
    slowest.X_end = max(slowest.X);
    slowest.T_diff_end = params.T_liq - slowest.T_end;
    slowest.X_diff_end = 1 - slowest.X_end;
    [slowest.z_unique, unique_inds_z] = unique_diff(slowest.z, params.chars.unique_thresh);%unique(slowest.z);
    slowest.T_unique = slowest.T(unique_inds_z); %NEED TO USE THIS FOR UNIQUE
    slowest.t_unique = slowest.t(unique_inds_z);
    [slowest.T_unique, unique_inds_T] = unique_diff(slowest.T_unique, params.chars.unique_thresh);
    slowest.z_unique = slowest.z_unique(unique_inds_T);
    slowest.z_unique(end) = params.reactor.H;
    slowest.t_end = slowest.t_unique(end);
    slowest.t_unique = slowest.t_unique(unique_inds_T);
    slowest.t_unique(end) = slowest.t_end;


    slowest.dTdts = (slowest.T_end - slowest.T)./(slowest.t(end) - slowest.t);
    slowest.dTdzs = (slowest.T_end - slowest.T)./(slowest.z(end) - slowest.z);

    %Iterate through each mass
    trail.Tis = cell(1, params.Nms);
    trail.zis = cell(1, params.Nms);
    trail.tis = cell(1, params.Nms);
    trail.uis = cell(1, params.Nms);

    %Update odeparams
    im = params.gparams.im;
    m = params.mms(im);
    n_i = m/params.M_gas_i;
    n_Ar_i = n_i * params.X_Ar;
    n_CH4_i = n_i * (1-params.X_Ar);
    odeparams.im = im;
    odeparams.m = m;
    %odeparams.Xi = 0; %
    odeparams.n_i = n_i;
    odeparams.n_Ar_i = n_Ar_i;
    odeparams.n_CH4_i = n_CH4_i;

    %Calculate slowest characteristics
    [T1,Y1] = TemperatureCharacteristic(0, params.t_f, T_i, [], params, odeparams, false);

    %Determine how many trailing points to have 
    t = T1;
    u = Y1(:,1);
    z = Y1(:,2);
    T = Y1(:,3);
    X = Y1(:,4);
    t_end_ind = min(find(z >=  params.reactor.H));
    t_res = T1(t_end_ind);
    T_end = max(T);
    X_end = max(X);
    T_diff = T_end - slowest.T_end;
    N_trail_pts = ceil(T_diff/params.dT_end_max);

    %Debug plot
    if params.chars.debug
        f1 = figure();

        % subplot(1,2,1);
        % plot(slowest.T, slowest.t); hold on;
        % plot(T, t);

        %subplot(1,2,2);
        p1 = plot(slowest.T, slowest.z, 'Linewidth', 1.5); hold on;
        plot(T, z, 'LineWidth', 1.5);
    end

    %Calculate initial temperature and 
    switch params.LC_trail_type
        case 'Final'
            [Tis, zis, chars] = TrailingPointsFinal(T1, Y1, slowest, params, odeparams);
        case 'Initial'
            [Tis, zis, tis, uis] = TrailingPointsInitial(T1, Y1, slowest, n_i,  params, odeparams);
    end

    %Apply 
    xlabel('Temperature (K)'); ylabel('z-position (m)');
    title('Trailing Characteristics')
    axis square; grid on; grid minor; 
    set(gca, 'FontSize', 18)

    %Log results.
    trail.Tis = fliplr(Tis);
    trail.zis = fliplr(zis);
    trail.tis = fliplr(tis);
    trail.uis = fliplr(uis);
    trail.slowest = slowest;

    %Close figure
    if params.chars.debug
        close(f1);
    end


        % %Calculate final temperatures - bias toward 
        % if N_trail_pts == 1
        %     Tf_trails = slowest.T_end; %The trail point is just the slowest case
        % else
        %     T_diff_nearest = T_diff - (N_trail_pts - 1) .* params.dT_end_max; %The final T difference from the closest core point
        %     Tf_trails = logspace(log10(slowest.T_end), log10(slowest.T_end + (N_trail_pts - 1) * params.dT_end_max), N_trail_pts); %slowest.T_end:params.dT_end_max:(slowest.T_end + (N_trail_pts - 1) * params.dT_end_max);
        % end
        

        % %Debug plot
        % if params.debug
        %     figure()
        % 
        %     subplot(1,2,1);
        %     plot(slowest.T, slowest.t); hold on;
        %     plot(T, t);
        % 
        %     subplot(1,2,2);
        %     plot(slowest.T, slowest.z); hold on;
        %     plot(T, z);
        % end
        % 
        % %Determine beginning z or t for trailing point engagement
        % dTdts = (T_end - T)./(t(end) - t);
        % dTdzs = (T_end - T)./(z(end) - z);
        % 
        % 
        % Tis = zeros(size(Tf_trails));
        % for ip = 2:length(Tf_trails)
        % 
        %     Tf_slow = slowest.T_end;
        %     Tf_b = Tf_trails(ip);
        % 
        %     %Plot point
        %     if params.debug
        %         subplot(1,2,2); 
        %         plot(Tf_b, params.reactor.H, 'ro'); 
        %     end
        % 
        %     %Identify Ti - starting point where two curves collide/diverge
        %     Ti_func_find = @(T_i) Ti_func(T_i, Tf_b, t_res, n_i, slowest, params, odeparams);

        %     frac = (1 - (Tf_b - slowest.T_end)/T_diff);
        %     Ti_guess = params.gparams.T_i + (slowest.T_end - T_i) .* (1 - (Tf_b - slowest.T_end)/T_diff);
        %     Tis_guess = zeros(1,1000);
        %     Tis_guess(1) = Ti_guess;
        %     err = 1;
        %     it = 1;
        %     err_tol = 0.1;
        %     while abs(err) > 0.1
        % 
        %         %Resolve 
        %         [Tf_act, err, T2, Y2] = Ti_func_find(Ti_guess);
        %         T_act = Y2(:,3);
        %         z_act = Y2(:,2);
        %         zi = z_act(1);
        % 
        %         %Plot
        %         if params.debug 
        %             p1 = plot(T_act, z_act);
        % 
        %         end
        % 
        %         %Update starting point
        %         Tis_guess_diffs = Ti_guess - Tis_guess(Tis_guess > 0);
        %         if err > 0
        %             negs = Tis_guess_diffs < 0;
        %             min_neg = min(abs(Tis_guess_diffs(negs)));
        %             if ~isempty(min_neg)
        %                 Ti_next_ind = find(Tis_guess_diffs == -min_neg);
        %                 Ti_next = Tis_guess(Ti_next_ind);
        % 
        % 
        %                 %[Ti_next, Ti_next_ind] = min(Tis_guess(Tis_guess_diffs < 0));
        %                 %Ti_next = Tis_guess(min(find(Tis_guess > Ti_guess)));
        %             %if ~isempty(Ti_next)
        %                 Ti_guess = (Ti_guess + Ti_next)/2;
        %             else
        %                 Ti_guess = Ti_guess + 0.5 * (Tf_slow - Ti_guess);
        %             end
        % 
        %         else
        %             poss = Tis_guess_diffs > 0;
        %             min_pos = min(abs(Tis_guess_diffs(poss)));
        %             if ~isempty(min_pos)
        %                 Ti_next_ind = find(Tis_guess_diffs == min_pos);
        %                 Ti_next = Tis_guess(Ti_next_ind);
        % 
        %                 %[~, Ti_next_ind] = min(Tis_guess_diffs(Tis_guess_diffs > 0));
        %                 %Ti_next = Tis_guess(Ti_next_ind)
        %                 %Ti_next = Tis_guess(max(find(Tis_guess < Ti_guess & Tis_guess > 0)));
        %                 %if ~isempty(Ti_next)
        %                 Ti_guess = (Ti_guess + Ti_next)/2;
        %             else
        %                 Ti_guess = Ti_guess - 0.5 * (Ti_guess - T_i);
        %             end
        % 
        %         end
        % 
        %         %Store guess
        %         it = it + 1;
        %         Tis_guess(it) = Ti_guess;
        % 
        %         %Plot result
        %         if params.debug & abs(err) > 0.1
        %             delete(p1);
        %         end
        % 
        %     end
        % 
        %     %Log starting points
        %     zis(ip) = zi;
        %     Tis(ip) = Ti_guess;
        % 
        % 
        % end


            
        %Check that it is correct
 

        % 
        % %Define coefficient 
        % 
        % 
        % p_bar = (params.p + params.p_surf)/2; %Mean pressure during rise
        % T_bar = (params.T_liq + params.T_i)/2; %Mean temperature during evolution
        % M_bar = params.M_gas; %Mean gas molar mass - Not accounting for reaciton
        % A_min_bar = params.A_min_func(T_bar, p_bar, M_bar);
        % A_min_theta = @(theta,p,M) params.A_min_func(params.T_liq - theta,p,M);
        % A_theta = @(theta,p,M) params.A_func(params.T_liq - theta,p,M);
        % 
        % %Define and solve ODE
        % theta_i = params.T_liq - params.T_i;
        % theta_largest_ode = @(t,theta) -A_min_theta(theta,p_bar,M_bar)*theta;
        % theta_ode = @(t,theta) -A_theta(theta,p_bar,M_bar) * theta;
        % ts_in = 0:params.dt:params.t_f;
        % [ts, thetas_largest] = ode45(theta_largest_ode, ts_in, theta_i); %For largest bubble
        % [ts, thetas] = ode45(theta_ode, ts_in, theta_i); %For reference bubble 
        % 
        % 
        % 
        % %Iterate through timesteps and determine the temperature difference
        % DTs = thetas_largest - thetas;
        % DT_max = max(DTs);


    

end



function [Tf_est, Tf_diff, T2, Y2] = Ti_func(T_i, Tf, t_res, n_i, slowest, params, odeparams)

    %Define maximum t based on residence time 
    t_f = t_res + 0.01;

    %Determine initial z

    zi = interp1(slowest.T_unique, slowest.z_unique, T_i);
    p_i = params.p_z(zi);

    %Estimate initial ui (assume terminal velocity)
    Vb = (n_i .* params.R .* T_i)/p_i;
    Db = (6 * Vb/pi)^(1/3);
    T_L = params.T_Lz(zi);
    rho_L = params.liquid.density(T_L);
    mu_L = params.liquid.dyn_visc(T_L);
    sigma = params.liquid.surf_tension(T_L);
    rho_G = 0.5;
    Eob = ((rho_L - rho_G) .* params.g .* Db.^2)./sigma;
    ui = CalcVelocity(Db, rho_L, mu_L, sigma, Eob, params.fsolve_opts);

    %Specify initial conversion
    Xi = odeparams.Xi;

    %Run simulation
    y0 = [ui; zi; T_i; Xi];
    [T2,Y2] = TemperatureCharacteristic(0, t_f, T_i, y0, params, odeparams, params.chars.debug);

    %Pull final values
    T = Y2(:,3);
    Tf_est = max(T);
    Tf_diff = Tf_est - Tf; %Shoould become zero

end






function [Tis, zis] = TrailingPointsFinal(T1, Y1, slowest, params, odeparams)

    %Determine how many trailing points to have 
    t = T1;
    u = Y1(:,1);
    z = Y1(:,2);
    T = Y1(:,3);
    X = Y1(:,4);
    t_end_ind = min(find(z >=  params.reactor.H));
    t_res = T1(t_end_ind);
    T_end = max(T);
    X_end = max(X);
    T_diff = T_end - slowest.T_end;
    N_trail_pts = ceil(T_diff/params.dT_end_max);


    %Calculate final temperatures - bias toward 
    if N_trail_pts == 1
        Tf_trails = slowest.T_end; %The trail point is just the slowest case
    else
        T_diff_nearest = T_diff - (N_trail_pts - 1) .* params.dT_end_max; %The final T difference from the closest core point
        Tf_trails = logspace(log10(slowest.T_end), log10(slowest.T_end + (N_trail_pts - 1) * params.dT_end_max), N_trail_pts); %slowest.T_end:params.dT_end_max:(slowest.T_end + (N_trail_pts - 1) * params.dT_end_max);
    end
    
    
    

    %Determine beginning z or t for trailing point engagement
    dTdts = (T_end - T)./(t(end) - t);
    dTdzs = (T_end - T)./(z(end) - z);


    Tis = zeros(size(Tf_trails));
    for ip = 2:length(Tf_trails)
        
        Tf_slow = slowest.T_end;
        Tf_b = Tf_trails(ip);

        %Plot point
        if params.debug & false
            subplot(1,2,2); 
            plot(Tf_b, params.reactor.H, 'ro'); 
        end

        %Identify Ti - starting point where two curves collide/diverge
        Ti_func_find = @(T_i) Ti_func(T_i, Tf_b, t_res, n_i, slowest, params, odeparams);

        frac = (1 - (Tf_b - slowest.T_end)/T_diff);
        Ti_guess = params.gparams.T_i + (slowest.T_end - T_i) .* (1 - (Tf_b - slowest.T_end)/T_diff);
        Tis_guess = zeros(1,1000);
        Tis_guess(1) = Ti_guess;
        err = 1;
        it = 1;
        err_tol = 0.1;
        while abs(err) > 0.1

            %Resolve 
            [Tf_act, err, T2, Y2] = Ti_func_find(Ti_guess);
            T_act = Y2(:,3);
            z_act = Y2(:,2);
            zi = z_act(1);

            %Plot
            if params.debug 
                p1 = plot(T_act, z_act);
                
            end

            %Update starting point
            Tis_guess_diffs = Ti_guess - Tis_guess(Tis_guess > 0);
            if err > 0
                negs = Tis_guess_diffs < 0;
                min_neg = min(abs(Tis_guess_diffs(negs)));
                if ~isempty(min_neg)
                    Ti_next_ind = find(Tis_guess_diffs == -min_neg);
                    Ti_next = Tis_guess(Ti_next_ind);
                

                    %[Ti_next, Ti_next_ind] = min(Tis_guess(Tis_guess_diffs < 0));
                    %Ti_next = Tis_guess(min(find(Tis_guess > Ti_guess)));
                %if ~isempty(Ti_next)
                    Ti_guess = (Ti_guess + Ti_next)/2;
                else
                    Ti_guess = Ti_guess + 0.5 * (Tf_slow - Ti_guess);
                end

            else
                poss = Tis_guess_diffs > 0;
                min_pos = min(abs(Tis_guess_diffs(poss)));
                if ~isempty(min_pos)
                    Ti_next_ind = find(Tis_guess_diffs == min_pos);
                    Ti_next = Tis_guess(Ti_next_ind);

                    %[~, Ti_next_ind] = min(Tis_guess_diffs(Tis_guess_diffs > 0));
                    %Ti_next = Tis_guess(Ti_next_ind)
                    %Ti_next = Tis_guess(max(find(Tis_guess < Ti_guess & Tis_guess > 0)));
                    %if ~isempty(Ti_next)
                    Ti_guess = (Ti_guess + Ti_next)/2;
                else
                    Ti_guess = Ti_guess - 0.5 * (Ti_guess - T_i);
                end

            end

            %Store guess
            it = it + 1;
            Tis_guess(it) = Ti_guess;

            %Plot result
            if params.debug & abs(err) > 0.1
                delete(p1);
            end

        end
        
        %Log starting points
        zis(ip) = zi;
        Tis(ip) = Ti_guess;


    end




end



function [Tis, zis, tis, uis] = TrailingPointsInitial(T1, Y1, slowest, n_i, params, odeparams)

    interp_method = params.chars.interp_method;

    %Pull values
    t = T1;
    u = Y1(:,1);
    z = Y1(:,2);
    T = Y1(:,3);
    X = Y1(:,4);

    %Determine how many trailing points to have  
    t_end_ind = min(find(z >=  params.reactor.H));
    if isempty(t_end_ind)
        t_end_ind = length(t);
    end

    t_res = T1(t_end_ind);
    T_end = max(T);
    X_end = max(X);
    
    slowest.DeltaT = slowest.T_end - slowest.T(1);
    T_diff = T_end - slowest.T_end;
    N_trail_pts = ceil(T_diff/params.dT_end_max);
    %N_trail_pts = ceil(slowest.DeltaT/params.dT_end_max);
    N_trail_pts = max([N_trail_pts, 2]);

    %Select 
    switch params.LC_trail_spacing
        case 'log'
            Tis = logspace(log10(params.T_min), log10(slowest.T_end), N_trail_pts+2);
            Tis = Tis(2:end-1);
        case 'lin'
            Tis = linspace(params.T_min, slowest.T_end, N_trail_pts+2);
            Tis = Tis(2:end-1);
        case 'geo'
            r_sum = params.LC_trail_geo_r.^[1:N_trail_pts+1] ;
            deltaT_i = (slowest.T_end - params.T_min)/sum(r_sum);
            deltaT_is = r_sum .* deltaT_i;
            deltaT_i_cumsum = cumsum(deltaT_is);
            Tis = [params.T_min, params.T_min + deltaT_i_cumsum];
            Tis = Tis(2:end-1);

    end

    %Find starting point
    zis = zeros(size(Tis));
    tis = zis;
    uis = zis;
    for ip = 1:N_trail_pts

        %Pull T_i
        Ti = Tis(ip);
    
        %Determine initial 
        zis(ip) = interp1(slowest.T_unique, slowest.z_unique, Ti);
        p_i = params.p_z(zis(ip));

        %Determine initial t
        tis(ip) = interp1(slowest.T_unique, slowest.t_unique, Ti);

        %Estimate initial ui (assume terminal velocity)
        Vb = (n_i .* params.R .* Ti)/p_i;
        Db = (6 * Vb/pi)^(1/3);
        T_L = params.T_Lz(zis(ip));
        rho_L = params.liquid.density(T_L);
        mu_L = params.liquid.dyn_visc(T_L);
        sigma = params.liquid.surf_tension(T_L);
        rho_G = 0.5;
        Eob = ((rho_L - rho_G) .* params.g .* Db.^2)./sigma;
        uis(ip) = CalcVelocity(Db, rho_L, mu_L, sigma, Eob, params.fsolve_opts);

    end


    %Resolve characteristics for each
    %Tis = zeros(1,N_trail_pts);
    if params.chars.debug & true
        %zis = Tis;
        for ip = 1:N_trail_pts
    
            %Pull T_i
            Ti = Tis(ip);
    
            %Define maximum t based on residence time 
            t_f = t_res + 0.01;
    
            %Determine initial z
            zi = interp1(slowest.T_unique, slowest.z_unique, Ti);
            p_i = params.p_z(zi);
    
            %Estimate initial ui (assume terminal velocity)
            Vb = (n_i .* params.R .* Ti)/p_i;
            Db = (6 * Vb/pi)^(1/3);
            T_L = params.T_Lz(zi);
            rho_L = params.liquid.density(T_L);
            mu_L = params.liquid.dyn_visc(T_L);
            sigma = params.liquid.surf_tension(T_L);
            rho_G = 0.5;
            Eob = ((rho_L - rho_G) .* params.g .* Db.^2)./sigma;
            ui = CalcVelocity(Db, rho_L, mu_L, sigma, Eob, params.fsolve_opts);
    
            %Specify initial conversion
            Xi = odeparams.Xi;
    
            %Run simulation
            y0 = [ui; zi; Ti; Xi];
            [T2,Y2] = TemperatureCharacteristic(0, t_f, Ti, y0, params, odeparams, false);
    
            %Pull final values
            z = Y2(:,2);
            T = Y2(:,3);
            Tfs(ip) = max(T);
    
            %Plot characteristic
            if params.debug & true
                plot(T, z, 'Linewidth', 1.5);
            end
    
            % %Store characteristics
            % zis(ip) = zi;
            % 
            % inds = 1:1:min(find(Y2(:,2) > params.reactor.H));
            % ts = T2(inds);
            % Zs = z(inds);
            % Ts = T(inds);
            % Xs = Y2(inds,4);
            % Zs_T{it}{ix} = Zs;
            % Ts_T{it}{ix} = Ts;
            % Xs_T{it}{ix} = Xs;
            % ts_T{it}{ix} = ts;
            % 
            % chars{im}.trail.Ts_zs{it}{ix} = interp1(Zs, Ts, zs_interp, interp_method);
            % chars{im}.trail.Ts_bs{it}{ix} = interp1(Zs, Ts, mesh.yy(:,1), interp_method);
            % chars{im}.trail.Ts_cs{it}{ix} = interp1(Zs, Ts, mesh.volcell_cents(:,2), interp_method);
            % chars{im}.trail.Xs_bs{it}{ix} = interp1(Zs, Xs, mesh.yy(:,1), interp_method);
            % chars{im}.trail.Xs_cs{it}{ix} = interp1(Zs, Xs, mesh.volcell_cents(:,2), interp_method);
            % chars{im}.trail.Xs_zs{it}{ix} = interp1(Zs, Xs, zs_interp, interp_method);
            % chars{im}.trail.ts_zs{it}{ix} = interp1(Zs, ts, zs_interp, interp_method);
    
    
    
            %Tis(ip) = Tis;
    
    
    
    
        end
    end

    

end

