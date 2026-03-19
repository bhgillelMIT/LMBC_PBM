%A function to resolve temperature and conversion characteristics for the
%PBM
function outputs = TemperatureCharacteristics(params)

    %
    if nargin < 1
        PBM_v3();
    end

    %Settings
    debug = params.chars.debug;
    print = params.chars.print;
    genplot = false;
    genplot_final = true;
    saveplot = true;
    lw = 1;
    ms = 12;
    fs = 18;
    dir_char = params.chars.fig_dir; % 'Figures/Characteristics/';
    dir_temp = [dir_char, 'T/'];
    dir_conv = [dir_char, 'X/']; %'Figures/Characteristics/X/';
    dir_data = params.chars.data_dir;  %'Data/Characteristics/';

    %Specify general plot settings
    plots.ms = ms; plots.lw = lw; plots.fs = fs;

    %Sepcify string corresponding to case
    outputs.filename = sprintf('Chars_%s_Nm=%d_Ti=%d_Tstd=%d_Tl=%d.mat', params.liquid.name, params.Nms, round(params.T_mu_i), params.T_std_i, round(params.T_liq));

    %Initial print
    if print
        fprintf('Characeristics - Beginning:\n')
    end

    %Make directories
    if saveplot

        %Create directory name
        dirname = 'Time_' + string(datetime);
        dirname = strrep(dirname, ' ', '_'); dirname = strrep(dirname, ':', '-');

        %Make time directory
        mkdir(dir_char + dirname);

        %Update directories
        dir_time = dir_char + dirname + '/';
        dir_temp = dir_char + dirname + '/T/';
        dir_conv = dir_char + dirname + '/X/';

        %Make T and X directories
        mkdir(dir_temp); mkdir(dir_conv);


    end

    %Colors
    colors.ChiliRed = [194/255, 24/255, 7/255];

    %Pull reference state
    T0 = params.T0;
    p0 = params.p0;

    %Pull materials
    methane = params.methane;
    hydrogen = params.hydrogen;
    carbon = params.carbon;

    %Pull mesh boundaries and 
    mesh = params.mesh;
    p_z = params.p_z; 

    zs_interp = 0:0.001:params.reactor.H;

    %Initialize ode parameter structure
    odeparams.H_func = @(n_CH4, n_H2, n_C, T) n_CH4 * integral(methane.Cp, T0, T) + ...
        n_H2 * integral(hydrogen.Cp, T0, T) + n_C * integral(carbon.Cp, T0, T);
    odeparams.p_z = p_z;
    odeparams.h_reactor = params.reactor.H;
    odeparams.T_Lz = params.T_Lz;
    odeparams.rho_L_T = params.rho_L_T;
    odeparams.g = params.g;
    odeparams.R = params.R;
    odeparams.methane = params.methane;
    odeparams.hydrogen = params.hydrogen;
    odeparams.carbon = params.carbon;
    odeparams.liquid = params.liquid;
    odeparams.argon = params.argon;
    odeparams.k0 = 6.6E13;
    odeparams.Ea = -370000;
    odeparams.h_conv = 5;
    odeparams.Xis = linspace(0,1, params.chars.N_Xis+1);%0:0.125:0.875;
    odeparams.Xis = odeparams.Xis(1:end-1);
    odeparams.h_conv_const = params.h_conv_const;
    odeparams.chars = params.chars;

    %Determine 
    m = params.mms(end); odeparams.m = m;
    n_i = m/params.M_gas_i;n_Ar_i = n_i * params.X_Ar; n_CH4_i = n_i * (1-params.X_Ar);
    odeparams.im = params.Nms;
    odeparams.Xi = 0;
    odeparams.n_i = n_i;
    odeparams.n_Ar_i = n_Ar_i;
    odeparams.n_CH4_i = n_CH4_i;
    T_initials = params.T_initials{end}; Tps = T_initials.Tps;
    [T,Y] = TemperatureCharacteristic(0, params.t_f, params.T_min, [], params, odeparams, true);
    Tf_min_largest = Y(end,3);

    %Iterate through each mass
    for im = 1:params.Nms

        %Pull mass
        m = params.mms(im);
        
        n_i = m/params.M_gas_i;
        n_Ar_i = n_i * params.X_Ar;
        n_CH4_i = n_i * (1-params.X_Ar);

        %Pull initial mass fraction
        init_fracs = params.T_initials{im}.Fracs;


        %Update ode parameters
        odeparams.im = im;
        odeparams.m = m;
        odeparams.Xi = 0; %
        odeparams.n_i = n_i;
        odeparams.n_Ar_i = n_Ar_i;
        odeparams.n_CH4_i = n_CH4_i;

        %Pull initial temperatures
        T_initials = params.T_initials{im};
        Tps = T_initials.Tps;
        NTs = length(Tps);

        icore = find(T_initials.engaged);
        itrail = find(~T_initials.engaged); 
        itrail = fliplr(itrail);

        %Allocate storage vectors
        Zs_T = cell(1, length(Tps));
        Ts_T = cell(1, length(Tps));
        Xs_T = cell(1, length(Tps));
        ts_T = cell(1, length(Tps));
        Ts_zs_T = cell(1, length(Tps));
        Ts_bs_T = cell(1, length(Tps));
        Ts_cs_T = cell(1, length(Tps));
        Xs_bs_T = cell(1, length(Tps));
        Xs_cs_T = cell(1, length(Tps));
        Xs_zs_T = cell(1, length(Tps));
        ts_zs_T = cell(1, length(Tps));

        %Calculate final temperature difference between this bubble size
        %and the largest bubble size
        [T,Y] = TemperatureCharacteristic(0, params.t_f, Tps(1), [], params, odeparams, debug);
        Tf_min_i = Y(end,3);
        DTf_i = Tf_min_i - Tf_min_largest;
        

        %Iterate through engaged pivots
        for it = icore 

            %Allocate storage array
            sol = struct();

            %Iterate through conversions
            for ix = 1:length(odeparams.Xis)

                %Set Xi
                odeparams.Xi = odeparams.Xis(ix);

                %Solve ODE system
                [T,Y] = TemperatureCharacteristic(0, params.t_f, Tps(it), [], params, odeparams, debug);
    
                %Interpolate results to the spatial grid
                end_ind = min(find(Y(:,2) > params.reactor.H));
                if isempty(end_ind)
                    end_ind = length(Y(:,2));
                end
                inds = 1:1:end_ind;
                ts = T(inds);
                Zs = Y(inds,2);
                Ts = Y(inds,3);
                Xs = Y(inds,4);


                %Store
                interp_method = params.chars.interp_method; %'linear';  
                Zs_T{it}{ix} = Zs;
                Ts_T{it}{ix} = Ts;
                Xs_T{it}{ix} = Xs;
                ts_T{it}{ix} = ts;
    
                Ts_zs_T{it}{ix} = interp1(Zs, Ts, zs_interp, interp_method);
                Ts_bs_T{it}{ix} = interp1(Zs, Ts, mesh.yy(:,1), interp_method);
                Ts_cs_T{it}{ix} = interp1(Zs, Ts, mesh.volcell_cents(:,2), interp_method);
                Xs_bs_T{it}{ix} = interp1(Zs, Xs, mesh.yy(:,1), interp_method);
                Xs_cs_T{it}{ix} = interp1(Zs, Xs, mesh.volcell_cents(:,2), interp_method);
                Xs_zs_T{it}{ix} = interp1(Zs, Xs, zs_interp, interp_method);
                ts_zs_T{it}{ix} = interp1(Zs, ts, zs_interp, interp_method);

                %Store initial fractions

                %Catch NaNs
                if any(isnan(Ts_cs_T{it}{ix}))
                    x = 1;
                end

            end

            

        end

        %Iterate through trailing points
        trail = T_initials.trail;
        for it = itrail

            %Determine the time when the point activates
            ui = trail.uis(it);
            Ti= trail.Tis(it);
            zi = trail.zis(it);
            t_active = interp1(trail.slowest.z_unique, trail.slowest.t_unique, zi); %s - moment when point activates
            if isnan(t_active)
                error('Invalid t_active');
            end
            
            %Calculate initial conversion
   
            [z_unique, unique_inds] = unique(trail.slowest.z);
            t_unique = trail.slowest.t(unique_inds);
            X_unique = trail.slowest.X(unique_inds);
            Xi = interp1(t_unique, X_unique, t_active);

            %Determine how many of the conversions to use
            Xis_below = odeparams.Xis < Xi;
            NXis_below = length(find(Xis_below));

            %Pull other unique values
            T_unique = trail.slowest.T_unique;
            u_unique = trail.slowest.u(unique_inds);   
            %t_unique = trail.slowest.t_unique;

            %Pull slowest characteristic off which the 
            N_pts = 200;
            ts_interp = linspace(0,t_active, N_pts);
            T1 = ts_interp';
            zs_intrp = interp1(t_unique, z_unique, T1);
            Ts_interp = interp1(trail.slowest.t_unique, T_unique, T1);
            us_interp = interp1(t_unique, u_unique, T1);
            Xs_interp = interp1(t_unique, X_unique, T1);
            Y1 = [us_interp, zs_intrp, Ts_interp, Xs_interp];

            %Iterate through conversions
            for ix = 1:length(odeparams.Xis)

                %Set Xi
                odeparams.Xi = max([Y1(end,end), odeparams.Xis(ix)]);
    
                %Log original odeparams
                odeparams_orig = odeparams;
    
                % T_i = Tps(it);
                % 
                % T_activate = T_i + T_initials.trail.dT;
                % 
                % %Determine the time when the point activates
                % try
                %     [Ts_interp, unique_inds] = unique(Ts_T{it+1}{ix}); %Filter out repetition because of interp function requirements
                %     ts_interp = ts_T{it+1}{ix};
                %     ts_interp = [ts_interp(unique_inds(1:end-1)); ts_interp(end)];
                %     t_active = interp1(Ts_interp, ts_interp, T_activate);
                % catch
                %     x = 1;
                % end

               
                

                % 
                % 
                % 
                % %Solve first part - at velocity of largest bubble
                % m_largest = params.mms(end);
                % n_i = m_largest/params.M_gas_i;
                % n_Ar_i = n_i * params.X_Ar;
                % n_CH4_i = n_i * (1-params.X_Ar);
                % 
                % odeparams.m = m_largest;
                % %odeparams.Xi = 0; %
                % odeparams.n_i = n_i;
                % odeparams.n_Ar_i = n_Ar_i;
                % odeparams.n_CH4_i = n_CH4_i;
                % [T1,Y1] = TemperatureCharacteristic(0, t_active, T_i, [], params, odeparams, debug);
    
                %Solve second part - at velocity of current bubble size
                T_i = Ti;
                y0 = [ui, zi, Ti, odeparams.Xi]; %Y1(end,:);
                odeparams = odeparams_orig;
                

                [T2,Y2] = TemperatureCharacteristic(t_active, params.t_f, T_i, y0, params, odeparams, debug);
                
                %Compile two results
                T = [T1(1:end-1); T2];
                Y = [Y1(1:end-1,:); Y2];
    
                %Log results
                end_ind = min(find(Y(:,2) > params.reactor.H));
                if isempty(end_ind)
                    end_ind = length(Y(:,2));
                end
                inds = 1:1:end_ind;
                ts = T(inds);
                Zs = Y(inds,2);
                Ts = Y(inds,3);
                Xs = Y(inds,4);
                interp_method = 'linear';  
                Zs_T{it}{ix} = Zs;
                Ts_T{it}{ix} = Ts;
                Xs_T{it}{ix} = Xs;
                ts_T{it}{ix} = ts;
    
                Ts_zs_T{it}{ix} = interp1(Zs, Ts, zs_interp, interp_method);
                Ts_bs_T{it}{ix} = interp1(Zs, Ts, mesh.yy(:,1), 'linear');
                Ts_cs_T{it}{ix} = interp1(Zs, Ts, mesh.volcell_cents(:,2), 'linear');
                Xs_bs_T{it}{ix} = interp1(Zs, Xs, mesh.yy(:,1), 'linear');
                Xs_cs_T{it}{ix} = interp1(Zs, Xs, mesh.volcell_cents(:,2), 'linear');
                Xs_zs_T{it}{ix} = interp1(Zs, Xs, zs_interp, interp_method);
                ts_zs_T{it}{ix} = interp1(Zs, ts, zs_interp, interp_method);

            end
            

        end

        %Log results
        itrails{im} = itrail;
        Ts_zs{im} = Ts_zs_T;
        Ts_bs{im} = Ts_bs_T;
        Ts_cs{im} = Ts_cs_T;
        Xs_bs{im} = Xs_bs_T;
        Xs_cs{im} = Xs_cs_T;
        Xs_zs{im} = Xs_zs_T;
        DTf(im) = DTf_i;
        fracs_in{im} = init_fracs;


        %Plot characteristic plots
        if genplot & true
            
            %Specify initial concentration index to plot
            ix = 1;

            

            %Create figure
            figure('units', 'normalized','outerposition', [0,0,1,1]);
            subplot(1,2,1);
            plots.fs = fs;
            PlotTempCharacteristics(gca, 'T', mesh, zs_interp, itrail, ix, Ts_zs_T, Ts_bs_T, 0, im, plots)

            %Plot conversion characteristic (assuming initial zero)
            subplot(1,2,2)
            PlotConvCharacteristics(gca, 'fixed_Xi', mesh, zs_interp, itrail, 1, Xs_zs_T, Xs_bs_T, 0, im, plots)

            %Save figure
            if saveplot
                filename = sprintf('Chars_Base_T_X_im=%d.png',im);
                filepath = dir_time + filename;
                saveas(gcf,filepath)
            end

            


            %Plot full domain - for each initial concentration
            figure('units', 'normalized','outerposition', [0,0,1,1]);
            N_Xs = length(odeparams.Xis);
            N_axes = 2 * ceil(N_Xs/2); %Must be an even number
            plots.fs = 8;
            for ix = 1:length(odeparams.Xis)
                subplot(2,4,ix)
                PlotConvCharacteristics(gca, 'fixed_Xi', mesh, zs_interp, itrail, ix, Xs_zs_T, Xs_bs_T, 0, im, plots)
            end

            %
            %Save figure if requested
            if saveplot
                filename = sprintf('Chars_Conv_Xis_im=%d.png', im);
                filepath = dir_conv + filename;
                saveas(gcf,filepath)
            end

            %Plot full domain - for each temperature
            figure('units', 'normalized','outerposition', [0,0,1,1]);
            N_Ts = length(Xs_zs_T);
            plots.fs = 6;
            for it = 1:N_Ts
                subplot(4,7,it);
                Ti = Ts_zs_T{it}{1}(1);
                PlotConvCharacteristics(gca, 'fixed_Ti', mesh, zs_interp, itrail, it, Xs_zs_T, Xs_bs_T,  Ti, im, plots)
            end

            %Save figure if requested
            if saveplot
                filename = sprintf('Chars_Conv_Tis_im=%d.png', im);
                filepath = dir_conv + filename;
                saveas(gcf,filepath)
            end


        end
   
        %Print update
        fprintf('Mass index = %i - Characteristics Calculated.\n', im);



        close all



    end

    


    %Export characteristics
    chars.ind_order = {'m', 'Xi', 'Ti'};
    chars.itrails = itrails;
    
    chars.Ts_zs = Ts_zs;
    chars.Ts_bs = Ts_bs;
    chars.Ts_cs = Ts_cs;
    chars.Xs_bs = Xs_bs;
    chars.Xs_cs = Xs_cs;
    chars.Xs_zs = Xs_zs;
    chars.fracs_in = fracs_in; %

    %Merge bins
    params.T_z.chars = chars;
    merged = MergeBins(params);
    chars.merged = merged;
    params.T_z.chars = chars;

    %Save output file
    dir_out = strcat(dir_data, outputs.filename);
    save(dir_out, 'chars');
    outputs.chars = chars;
    
    %Plot full domain - individual
    if params.chars.debug
        figtemp = figure('units', 'normalized','outerposition', [0,0,1,1]);
        figconv = figure('units', 'normalized','outerposition', [0,0,1,1]);
        N_ms = length(Ts_zs);
        N_sps = 12; %Number of subplots
        N_figs = ceil(N_ms/N_sps);
        ix = 1;
        plots.fs = 8;
        for im = 1:N_ms
    
            rem = mod(im, N_sps);
    
            %Specify figure and subplot
            figure(figtemp)
            if rem == 0
                subplot(3,4,12);
            else
                subplot(3,4, rem);
            end
    
            %Pull initial value of conversion
            Xi = Xs_zs{im}{1}{1}(1);
    
            %Plot temperature characteristics
            Ts_zs_T = Ts_zs{im};
            Ts_bs_T = Ts_bs{im};
            itrail = itrails{im};
            PlotTempCharacteristics(gca, 'T', mesh, zs_interp, itrail, ix, Ts_zs_T, Ts_bs_T, Xi, im, plots)
    
            %Specify figure and subplot
            figure(figconv)
            if rem == 0
                subplot(3,4,12);
            else
                subplot(3,4, rem);
            end
    
            %Plot conversion characteristics
            Xs_zs_T = Xs_zs{im};
            Xs_bs_T = Xs_bs{im};
            PlotConvCharacteristics(gca, 'fixed_Xi', mesh, zs_interp, itrail, 1, Xs_zs_T, Xs_bs_T, 0, im, plots)
    
    
            %Create new figure
            if rem == 0
    
                %Save figure
                if saveplot
                    %Save temperature characteristics
                    filename = sprintf('Chars_T_ms_%d.png', ceil(im/N_sps));
                    filepath = dir_time + filename;
                    saveas(figtemp,filepath)
    
                    %Save conversion characteristics
                    filename = sprintf('Chars_X_ms_%d.png', ceil(im/N_sps));
                    filepath = dir_time + filename;
                    saveas(figconv,filepath)
    
                end
    
                %Create new figures
                figtemp = figure('units', 'normalized','outerposition', [0,0,1,1]);
                figconv = figure('units', 'normalized','outerposition', [0,0,1,1]);
    
            elseif im == N_ms
                if saveplot
                    %Save temperature characteristics
                    filename = sprintf('Chars_T_ms_%d.png', ceil(im/N_sps));
                    filepath = dir_time + filename;
                    saveas(figtemp,filepath)
    
                    %Save conversion characteristics
                    filename = sprintf('Chars_X_ms_%d.png', ceil(im/N_sps));
                    filepath = dir_time + filename;
                    saveas(figconv,filepath)
                end
            end
    
    
    
        end
    end

    %Plot conversion characteristics


%% Merge Bins

    

    


%% 3D temperature plot

    % %Create fixed mesh to use for surface
    % T_min = params.T_min; T_max = params.T_liq;
    % N_pts = 300;
    % Ts_grid = linspace(T_min, T_max, N_pts);
    % ms_grid = params.mms;
    % 
    % %Create plot
    % figure('units', 'normalized','outerposition', [0,0,1,1]);
    % 
    % ix = 1;
    % Ts_3D = []; ms_3D = []; zs_3D = [];
    % zs_z = mesh.yy(:,1);
    % for im = 1:N_ms
    % 
    %     %Pull mass
    %     m = params.mms(im);
    % 
    %     %Iterate through temperatures 
    %     N_Ts = length(Ts_bs{im});
    %     for it = 1:N_Ts
    %         Ts_line = Ts_bs{im}{it}{ix};
    %         ms_line = m .* ones(length(Ts_bs{im}{it}{ix}), 1);
    %         zs_line = zs_z; %repmat(zs_z, N_Ts, 1)
    % 
    % 
    %         %Plot 
    %         if im == 1 | im == N_ms
    %             plot3(Ts_line, ms_line, zs_line, 'r-'); hold on;
    %         end
    % 
    %         %Log first result
    %         if it == 1
    % 
    %             %Add points to represent top surface 
    %             if max(Ts_line) < T_max
    %                 Ts_new = linspace(max(Ts_line), T_max, 20);
    %                 zs_new = params.reactor.H .* ones(size(Ts_new'));
    %                 ms_new = m * ones(size(zs_new));
    %                 Ts_line = [Ts_line; Ts_new(2:end)'];
    %                 zs_line = [zs_line; zs_new(2:end)];
    %                 ms_line = [ms_line; ms_new(2:end)];
    %             end
    % 
    %             %Interpolate to a fixed temperatre grid
    %             try
    %                 zs_line = interp1(Ts_line, zs_line, Ts_grid);
    %             catch
    %                 x = 1;
    %             end
    %             Ts_line = Ts_grid;
    %             ms_line = m * ones(size(Ts_line));
    % 
    %             %Log values for
    %             ms_3D = [ms_3D; ms_line];
    %             zs_3D = [zs_3D; zs_line];
    %             Ts_3D = [Ts_3D; Ts_line];
    % 
    %         end
    % 
    %     end
    % 
    %     %Specify corresponding mass and z positions
    % 
    % end
    % 
    % 
    % %Create initial meshgrid and higher resolution grid
    % ms_grid_hi = linspace(min(params.mms), max(params.mms), 100);
    % [Ts_mesh, ms_mesh] = meshgrid(Ts_grid, ms_grid);
    % [Ts_mesh_hi, ms_mesh_hi] = meshgrid(Ts_grid, ms_grid_hi);
    % 
    % %Increase resolution
    % zs_mesh_hi = interp2(Ts_mesh, ms_mesh, zs_3D, Ts_mesh_hi, ms_mesh_hi, 'linear');
    % 
    % %Create surface plot
    % colormap autumn
    % surf(Ts_mesh_hi, ms_mesh_hi, zs_mesh_hi);
    % colorbar
    % 
    % %Plot aesthetics
    % m_min = min(params.mms); m_max = max(params.mms);
    % xlim([T_min, T_max]); ylim([m_min, m_max]); zlim([0, params.h])
    % xlabel('Temperature (C)'); ylabel('Mass (kg)'); zlabel('Z-Pos. (m)');
    
%% Create interpolation scheme 

    %
    


end


% function [T,Y] = TemperatureCharacteristic(t_i, t_f, Ti, y0, params, odeparams, debug)
% 
%     %
%     im = odeparams.im;
%     Xi = odeparams.Xi;
%     zi = 0;
%     ui = 0;
% 
%     %Define initial conditions
%     if isempty(y0)
%         y0 = [ui; zi; Ti; Xi];
%     end
% 
%     %Solve ODE system
%     t_i = round(t_i, 4);
%     options = odeset(RelTol=1e-8,AbsTol=1e-10);
%     tspan = [t_i, t_f]; %t_i:0.0001:t_f;
%     odefunc = @(t,y) TempOde(t,y,odeparams);
%     [T,Y] = ode45(odefunc, tspan, y0, options);
% 
%      %Debug plots
%     if debug
% 
%         %Plot values versus time
%         figure();
% 
%         subplot(2,2,1);
%         plot(T, Y(:,1), 'r-');
%         xlabel('Time (s)'); ylabel('Velocity (m/s)');
% 
%         subplot(2,2,2);
%         plot(T, Y(:,2), 'k-');
%         xlabel('Time (s)'); ylabel('Vertical Pos.');
% 
%         subplot(2,2,3);
%         plot(T, Y(:,3), 'r-');
%         xlabel('Time (s)'); ylabel('Temperature');
% 
%         subplot(2,2,4);
%         plot(T, Y(:,4), 'k-');
%         xlabel('Time (s)'); ylabel('Conversion');
% 
% 
%     end
% 
%     %Print update
%     fprintf('--Characteristic (i_m = %d; T_i = %0.2f K; X_i = %0.3f %%) - Done\n', im, Ti, 100.*Xi);
% 
% 
% 
% end
% 
% function dydt = TempOde(t,y,params)
% 
% 
% 
%     %Pull differential variables
%     u = y(1);
%     z = y(2);
%     T = y(3);
%     X = y(4);
% 
%     %Calculate moles from X
%     n_Ar = params.n_Ar_i;
%     n_CH4 = params.n_CH4_i * (1-X);
%     n_H2 = 2 * params.n_CH4_i * X;
%     n_C = 1 * params.n_CH4_i * X;
%     n_gas = n_CH4 + n_H2;
%     n_tot = n_gas + n_C;
% 
%     %Determine local conditions
%     T_L = params.T_Lz(z);
%     p = params.p_z(z);
%     V = (n_gas * params.R * T)/p;
%     Db = 2 * ((3 * n_gas * params.R * T)/(4 * pi * p))^(1/3);
%     As = 4 * pi * (Db/2).^2;
%     Ac = As/4; %m2 - cross sectional area
%     rho_L = params.rho_L_T(T_L);
%     rho_G = (p .* params.methane.mol_mass)./(params.R .* T);
%     mu_L = params.tin.dyn_visc(T);
%     sigma = params.tin.surf_tension(T);
%     Eob = ((rho_L - rho_G) .* params.g .* Db.^2)./sigma;
% 
%     %Calculate heat capacity
%     C = n_CH4 * params.methane.Cp(T) + n_H2 * params.hydrogen.Cp(T)...
%         + n_C * params.carbon.Cp(T) + n_Ar * params.argon.Cp(T);
% 
% 
%     %Calculate bubble velocity 
%     %ub = CalcVelocity(Db, rho_L, mu_L, sigma, Eob);
% 
%     %Calculate reaction rate
%     kr = params.k0 .* exp(params.Ea./(params.R .* T));
% 
%     %Forces on bubble
%     Fb = rho_L * V * params.g;
%     Fg = params.m * params.g;
%     Cd = CalcCd(u, Db, rho_L, mu_L, sigma, Eob);
%     Fd = 0.5 * rho_L * u^2 * Cd * Ac;
%     m_vm = (11/16) * rho_L * V;
% 
%     %Derivatives
%     if z < params.h_reactor
%         dudt = (Fb - Fd - Fg)/(m_vm + params.m);
%         dTdt = (params.h_conv * As)/(C) * (T_L - T);
%         dXdt = kr .* (1 - X);
%         dzdt = u;
%     else
%         dudt = 0; dTdt = 0; dXdt = 0; dzdt = 0;
%     end
% 
%     %Catch nans
%     if isnan(dudt)
%         x = 1;
%     end
% 
% 
%     %Define output vector
%     dydt = [dudt; dzdt; dTdt; dXdt];
% 
% 
% end
% 


function PlotTempCharacteristics(ax, mode, mesh, zs_interp, itrail, ix, Ts_zs_T, Ts_bs_T, Ti, im, plots)

    %Pull plot settings
    lw = plots.lw;
    ms = plots.ms;
    fs = plots.fs;

    %Plot horizontal lines to indicate spatial cells
    for iz = 1:length(mesh.yy(:,1))
        plot([800, 1300], [mesh.yy(iz,1), mesh.yy(iz,1)], 'k-', 'Color', [0.8, 0.8, 0.8]); hold on;
    end

    %Plot temperature lines
    T_floor = floor(min(Ts_zs_T{1}{1}-273.15)/100) * 100;
    T_ceil = ceil(max(Ts_zs_T{end}{end} - 273.15)/100) * 100;
    Ts_vert = T_floor:50:T_ceil; %[850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1250];
    for it = 1:length(Ts_vert)
        plot([Ts_vert(it), Ts_vert(it)], [0, 0.3],  'k-', 'Color', [0.8, 0.8, 0.8]); hold on;
    end

    %Iterate through temperatures
    for it = 1:length(Ts_bs_T)

        %Pull relevant values
        Ts = Ts_zs_T{it}{ix};
        zs = zs_interp;
        Ts_z = Ts_bs_T{it}{ix};
        zs_z = mesh.yy(:,1);
     
        %Plot trailing points
        if any(it == itrail)

            %Plot characteristic lines
            plot(Ts-273.15, zs, 'r-', 'LineWidth', lw, 'Color', [0.4, 0.4, 0.4]); hold on;

            %Plot points on spatial mesh
            
            plot(Ts_z-273.15, zs_z, 'r.', 'MarkerSize', ms, 'Color', [0.4, 0.4, 0.4]);

        else %Plot core and lead points

            %Plot characteristic lines
            plot(Ts-273.15, zs, 'r-', 'LineWidth', lw); hold on;

            %Plot points on spatial mesh  
            plot(Ts_z-273.15, zs_z, 'r.', 'MarkerSize', ms);

        end

        

    end

    %Aesthetics
    axis square;
    xlabel('Temperature (C)'); ylabel('Z-Position (m)');
    title(sprintf('Temp. Characteristics (X_i = %0.3f, i_m = %d)', Ti, im));
    set(gca, 'FontSize', fs, 'FontWeight', 'bold')
    %grid on; grid minor; 

    %Legend
    h = zeros(1,3);
    h(1) = plot(NaN, NaN, 'r-', 'LineWidth', lw);
    h(2) = plot(NaN, NaN, 'r-', 'LineWidth', lw, 'Color', [0.4, 0.4, 0.4]);
    h(3) = plot(NaN, NaN, 'r.', 'MarkerSize', ms, 'Color', [0.2, 0.2, 0.2]);
    legend(h, 'Core + Lead points', 'Trail points', 'Value at spatial boundaries', 'location', 'northwest', 'BackgroundAlpha', 0.7);
    


end


function PlotConvCharacteristics(ax, mode, mesh, zs_interp, itrail, i, Xs_zs_T, Xs_bs_T, Ti, im, plots)
    
    %Pull plot settings
    lw = plots.lw;
    ms = plots.ms;
    fs = plots.fs;

    %Pull initial Xi

    %Plot horizontal lines to indicate 
    for iz = 1:length(mesh.yy(:,1))
        plot([0,1], [mesh.yy(iz,1), mesh.yy(iz,1)], 'k-', 'Color', [0.8, 0.8, 0.8]); hold on;
    end

    %Plot temperature lines
    Xs_vert = 0.1:0.1:0.9;
    for it = 1:length(Xs_vert)
        plot([Xs_vert(it), Xs_vert(it)], [0, 0.3],  'k-', 'Color', [0.8, 0.8, 0.8]); hold on;
    end

    %
    if strcmp(mode, 'fixed_Xi')

        for it = 1:length(Xs_bs_T)
    
            %Pull relevant values
            Xs = Xs_zs_T{it}{i};
            zs = zs_interp;
            Xs_z = Xs_bs_T{it}{i};
            zs_z = mesh.yy(:,1);
    
            %Plot trailing points
            if any(it == itrail)
    
                %Plot characteristic lines
                plot(Xs, zs, 'r-', 'LineWidth', lw, 'Color', [0.4, 0.4, 0.4]); hold on;
    
                %Plot points on spatial mesh
                
                plot(Xs_z, zs_z, 'r.', 'MarkerSize', ms, 'Color', [0.4, 0.4, 0.4]);
    
            else %Plot core and lead points
    
                %Plot characteristic lines
                plot(Xs, zs, 'b-', 'LineWidth', lw); hold on;
    
                %Plot points on spatial mesh  
                plot(Xs_z, zs_z, 'b.', 'MarkerSize', ms);
    
            end
    
        end

        Xi = Xs_zs_T{it}{i}(1);

        %Aesthetics
        axis square;
        xlabel('Conversion (-)'); ylabel('Z-Position (m)');
        titlestr = sprintf('Conversion Characteristics (X_i = %0.3f; i_m = %i)', Xi, im);
        title(titlestr)
        set(gca, 'FontSize', fs, 'FontWeight', 'bold')
    
        %Legend
        h = zeros(1,3);
        h(1) = plot(NaN, NaN, 'r-', 'LineWidth', lw);
        h(2) = plot(NaN, NaN, 'r-', 'LineWidth', lw, 'Color', [0.4, 0.4, 0.4]);
        h(3) = plot(NaN, NaN, 'r.', 'MarkerSize', ms, 'Color', [0.2, 0.2, 0.2]);
    
        legend(h, 'Core + Lead points', 'Trail points', 'Value at spatial boundaries', 'location', 'northwest', 'FontSize', fs, 'BackgroundAlpha', 0.7);

    elseif strcmp(mode, 'fixed_Ti')

        for ix = 1:length(Xs_zs_T{i})

            %Pull relevant values
            Xs = Xs_zs_T{i}{ix};
            zs = zs_interp;
            Xs_z = Xs_bs_T{i}{ix};
            zs_z = mesh.yy(:,1);
    
            %Plot trailing points
            if any(i == itrail)
    
                %Plot characteristic lines
                plot(Xs, zs, 'r-', 'LineWidth', lw, 'Color', [0.4, 0.4, 0.4]); hold on;
    
                %Plot points on spatial mesh
                
                plot(Xs_z, zs_z, 'r.', 'MarkerSize', ms, 'Color', [0.4, 0.4, 0.4]);
    
            else %Plot core and lead points
    
                %Plot characteristic lines
                plot(Xs, zs, 'b-', 'LineWidth', lw); hold on;
    
                %Plot points on spatial mesh  
                plot(Xs_z, zs_z, 'b.', 'MarkerSize', ms);
    
            end

        end


        %Aesthetics
        axis square;
        xlabel('Conversion (-)'); ylabel('Z-Position (m)');
        titlestr = sprintf('Conv. Characteristics (T_i = %0.3f; i_m = %i)', Ti, im);
        title(titlestr)
        set(gca, 'FontSize', fs, 'FontWeight', 'bold')
    
        %Legend
        h = zeros(1,3);
        h(1) = plot(NaN, NaN, 'r-', 'LineWidth', lw);
        h(2) = plot(NaN, NaN, 'r-', 'LineWidth', lw, 'Color', [0.4, 0.4, 0.4]);
        h(3) = plot(NaN, NaN, 'r.', 'MarkerSize', ms, 'Color', [0.2, 0.2, 0.2]);
    
        legend(h, 'Core + Lead points', 'Trail points', 'Value at cell bounds.', 'location', 'northwest', 'FontSize', fs);



    end
end