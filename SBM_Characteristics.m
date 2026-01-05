function chars = SBM_Characteristics(params)

%% Setup

    %Debug
    debug = true;

    %Clean up
    close all

    %Plot settings
    lw = 2;
    fs = 18;
    ms = 28;

    %Load materials
    addpath('Material Properties/');
    load tin.mat
    load hydrogen.mat
    load carbon.mat
    load methane.mat
    load argon.mat

%% Load data

    if nargin < 1
        PBM_v3();
    end

    %Settings
    interp_method = 'linear'; 
    calc_T_MC = false;

    %Reference states
    T0 = 298.15;
    p0 = 101325;

    %Inputs
    folder = params.SBM_folder; %'SBM Characteristics/Demo2/';
    
    

    %Load files
    files = dir(folder);
    N_files = length(find([files.bytes] > 0));
    files = files(3:end);

    %Reorder by date
    nums = zeros(N_files, 1);
    for i = 1:N_files
        name = files(i).name;
        nums(i) = sscanf(name, 'Char_%d');
    end
    dates = files(:).date;
    [sorted_nums, sorted_inds] = sort(nums);
    files = files(sorted_inds);

    %Pull face and center locations from mesh
    mesh_bounds = params.mesh.yy(:,1);
    mesh_cents = params.mesh.volcell_cents(:,2);
    zs_interp = linspace(0, params.reactor.H, 2000);


    %Allocate storage arrays
    Ts_zs = zeros(N_files, length(zs_interp));
    Ts_bs = zeros(N_files, length(mesh_bounds)); 
    Ts_cs = zeros(N_files, length(mesh_cents)); 
    Xs_bs = zeros(N_files, length(mesh_bounds)); 
    Xs_cs = zeros(N_files, length(mesh_cents)); 
    Xs_zs = zeros(N_files, length(zs_interp));
    ts_zs = zeros(N_files, length(zs_interp));




    
    %Iterate through files
    N_files = length(files); %REMOVE REMOVE REMOVE REMOVE REMOVE
    Tis = zeros(N_files, 1);
    Xis = zeros(N_files, 1);
    ris = zeros(N_files, 1);
    mis = zeros(N_files, 1);
    T_MCs = cell(N_files, 1);
    T_bars = cell(N_files, 1);
    Xs = cell(N_files, 1);

    for i = 1:(N_files)

        %Load file
        file = files(i);
        filename = file.name;
        filepath = [folder, filename];
        data = readtable(filepath);

        %Isolate variables
        cnames = data.Properties.VariableNames;
        t = data.('TIME');
        z = data.('EXPERIMENTAL_B_H');
        h_reactor = z(end);
        ri = data.('EXPERIMENTAL_R_INITIAL');
        mi = data.('EXPERIMENTAL_M_INITIAL');
        Nt = length(t);
        T_bar = data.('EXPERIMENTAL_B_T_BAR');
        X = data.('EXPERIMENTAL_B_X_CONV');

        %Truncate based 

        %Determine initial conditions
        Xis(i) = X(1);
        Tis(i) = T_bar(1);
        mis(i) = mi(1);
        ris(i) = ri(1);

        %Pull all temperatures
        N_shells = length(find(contains(cnames, 'EXPERIMENTAL_B_TS_')));
        Ts = zeros(Nt, N_shells);
        Cps = Ts;
        ns = Ts;
        ns_CH4 = Ts; ns_Ar = Ts; ns_C = Ts; ns_H2 = Ts;
        for is = 1:N_shells
            Cpstr = sprintf('EXPERIMENTAL_B_CPS_%d_', is);
            Tstr = sprintf('EXPERIMENTAL_B_TS_%d_', is);
            ns_str = sprintf('EXPERIMENTAL_B_NS_%d_', is);
            ns_CH4_str = sprintf('EXPERIMENTAL_B_NS_CH4_%d_', is);
            ns_Ar_str = sprintf('EXPERIMENTAL_B_NS_AR_%d_', is);
            ns_H2_str = sprintf('EXPERIMENTAL_B_NS_H2_%d_', is);
            ns_C_str = sprintf('EXPERIMENTAL_B_NS_C_%d_', is);
            Ts(:,is) = data.(Tstr);
            Cps(:,is) = data.(Cpstr);
            ns(:,is) = data.(ns_str);
            ns_CH4(:,is) = data.(ns_CH4_str);
            ns_Ar(:,is) = data.(ns_Ar_str);
            ns_H2(:,is) = data.(ns_H2_str);
            ns_C(:,is) = data.(ns_C_str);
        end

        %Calculate mixing cup temperature
        if calc_T_MC
            T_MC = zeros(length(t), 1);
            for it = 1:length(t)
    
                %Calculate total number of moles in the bubble at time t(it)
                n_CH4 = sum(ns_CH4(it, :));
                n_Ar = sum(ns_Ar(it, :));
                n_H2 = sum(ns_H2(it, :));
                n_C = sum(ns_C(it, :));
                
                
                Hs = zeros(1,N_shells);
                for is = 1:N_shells
                    H_CH4 = ns_CH4(it,is) .* integral(methane.Cp, T0, Ts(it, is));
                    H_AR = ns_Ar(it,is) .* integral(argon.Cp, T0, Ts(it, is));
                    H_H2 = ns_H2(it,is) .* integral(hydrogen.Cp, T0, Ts(it, is)); 
                    H_C = ns_C(it,is) .* integral(carbon.Cp, T0, Ts(it, is)); 
                    Hs(is) = H_CH4 + H_H2 + H_AR + H_C;
                end
    
                %Calculate total enthalpy of bubble
                H_total = sum(Hs);
    
                %Calculate mixing cup temperature which represents average
                T_func = @(T) n_CH4 * integral(methane.Cp, T0, T) + n_Ar * integral(argon.Cp, T0, T)...
                    + n_H2 * integral(hydrogen.Cp, T0, T) + n_C * integral(carbon.Cp, T0, T) - H_total;
                T_MC(it) = fzero(T_func, T_bar(it));
    
            end

        else
            T_MC = T_bar;

        end

        %Calculate values at mesh points
        [z_uni, uni_inds] = unique(z);
        z_uni(1) = 0; z_uni(end) = params.reactor.H;
        Ts_zs(i,:) = interp1(z_uni, T_MC(uni_inds), zs_interp, interp_method);
        Ts_bs(i,:) = interp1(z_uni, T_MC(uni_inds), mesh_bounds, interp_method);
        Ts_cs(i,:) = interp1(z_uni, T_MC(uni_inds), mesh_cents, interp_method);
        Xs_bs(i,:) = interp1(z_uni, X(uni_inds), mesh_bounds, interp_method);
        Xs_cs(i,:) = interp1(z_uni, X(uni_inds), mesh_cents, interp_method);
        Xs_zs(i,:) = interp1(z_uni, X(uni_inds), zs_interp, interp_method);
        ts_zs(i,:) = interp1(z_uni, t(uni_inds), zs_interp, interp_method);

        %Store outpouts
        T_bars{i} = T_bar;
        T_MCs{i} = T_MC;
        Xs{i} = X;






    end


    %Identify unique starting points and sort them
    ms_unique = unique(mis); %Needs to be imposed 
    ris_unique = unique(ris); 
    Xis_unique = unique(Xis); Xis_unique = sort(Xis_unique);
    Tis_unique = unique(Tis); Tis_unique = sort(Tis_unique);
    N_Tis = length(Tis_unique);
    N_Xis = length(Xis_unique);
    N_ms = length(ris_unique);

    %Create characteristic structure
    chars = struct;
    chars.ms = cell(1, N_ms);
    %chars.Ts_zs = cell(1, N_Tis); chars.Ts_bs = cell(1, N_Tis); chars.Ts_cs = cell(1, N_Tis);
    %chars.Xs_bs = cell(1, N_Tis); chars.Xs_cs = cell(1, N_Tis); chars.Xs_zs = cell(1, N_Tis);

    %Create cell arrays
    Ts_cell = cell(1, N_Tis);
    Xs_cell = cell(1, N_Xis);
    for im = 1:N_ms

        %List initial temperatures
        chars.ms{im}.Tis = Tis_unique ; %Different initial temperatures for each masss
        chars.ms{im}.Xis = Xis_unique;

        %
        chars.ms{im}.Ts_zs = cell(1, N_Tis); chars.ms{im}.Ts_bs = cell(1, N_Tis); chars.ms{im}.Ts_cs = cell(1, N_Tis);
        chars.ms{im}.Xs_bs = cell(1, N_Tis); chars.ms{im}.Xs_cs = cell(1, N_Tis); chars.ms{im}.Xs_zs = cell(1, N_Tis);

        for it = 1:N_Tis
            chars.ms{im}.Ts_zs{it} = Xs_cell;
            chars.ms{im}.Ts_bs{it} = Xs_cell;
            chars.ms{im}.Ts_cs{it} = Xs_cell;
            chars.ms{im}.Xs_zs{it} = Xs_cell;
            chars.ms{im}.Xs_bs{it} = Xs_cell;
            chars.ms{im}.Xs_cs{it} = Xs_cell;
        end
    end
    

    %Sort files into final structure
    for i = 1:N_files

        %Determine index
        ri = ris(i);
        m = mis(i);
        Ti = Tis(i);
        Xi = Xis(i);

        %Determine temperatures and conversions at cell boundaries and
        %   centers
        r_ind = find(ris_unique == ri);
        m_ind = find(ms_unique == m);
        T_ind = find(Tis_unique == Ti);
        X_ind = find(Xis_unique == Xi);

        %Store in the characteristics array
        chars.ms{r_ind}.Ts_zs{T_ind}{X_ind} = Ts_zs(i,:);
        chars.ms{r_ind}.Ts_bs{T_ind}{X_ind} = Ts_bs(i,:);
        chars.ms{r_ind}.Ts_cs{T_ind}{X_ind} = Ts_cs(i,:);
        chars.ms{r_ind}.Xs_zs{T_ind}{X_ind} = Xs_zs(i,:);
        chars.ms{r_ind}.Xs_bs{T_ind}{X_ind} = Xs_bs(i,:);
        chars.ms{r_ind}.Xs_cs{T_ind}{X_ind} = Xs_cs(i,:);


    end


%% Plot results

    figure();

    for it = 1:N_Tis
        ix = 1
        Ts_z = chars.ms{1}.Ts_zs{it}{ix}
        Xs_z = chars.ms{1}.Xs_zs{it}{ix}

        subplot(1,2,1);
        plot(Ts_z, zs_interp, 'r-'); hold on;
        
        
        subplot(1,2,2);
        plot(Xs_z, 100.*zs_interp, 'k-'); hold on;

    end

    subplot(1,2,1);
    xlabel('Temperature (K)'); ylabel('Z-Position (cm)');
    grid on; grid minor; axis square;
    title('Temperature Curves');
    set(gca, 'FontSize', 18);

    subplot(1,2,2);
    xlabel('Temperature (K)'); ylabel('Z-Position (cm)');
    grid on; grid minor; axis square;
    title('Conversion Curves');
    set(gca, 'FontSize', 18);

    



%% Process data



    params = params;

end




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
    Ts_vert = [850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1250];
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

    %Legend
    h = zeros(1,3);
    h(1) = plot(NaN, NaN, 'r-', 'LineWidth', lw);
    h(2) = plot(NaN, NaN, 'r-', 'LineWidth', lw, 'Color', [0.4, 0.4, 0.4]);
    h(3) = plot(NaN, NaN, 'r.', 'MarkerSize', ms, 'Color', [0.2, 0.2, 0.2]);
    legend(h, 'Core + Lead points', 'Trail points', 'Value at spatial boundaries', 'location', 'northwest', 'BackgroundAlpha', 0.7);



end