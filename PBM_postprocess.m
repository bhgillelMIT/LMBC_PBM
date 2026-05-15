function results = PBM_postprocess(t, y, params)
% A function to process the outputs from the PBM. Analyzes the bubble size
% distribution, volumetric and mass flow rates, and other physical aspects
% of the simulation.


%% Setup

    %Handle figures
    %set(0, 'DefaultFigureWindowStyle', 'docked');

    %Handle no inputs
    if nargin < 1
        close all 
        addpath('src/')
        addpath('Data/Solutions/');
        output = load('PBM_T=1_X=1_10-May-2026_22-50-27.mat'); %load('PBM_T=1_X=0_23-Apr-2026_12-00-11.mat'); %load('Data/Solutions/PBM_output_20-Feb-2026_16-00-09.mat');
        output = output.output;
        t = output.T;
        y = output.Y;
        params = output.params;
        params.sol.postprocess = true;
    end

    %Check if postprocess boolean exists
    try
        params.sol.postprocess;
    catch
        params.sol.postprocess = true;
    end

    %Plot settings
    lw = 1.5;
    fs = 18;

    %Pull materials
    methane = params.methane;
    hydrogen = params.hydrogen;
    argon = params.argon;
    carbon = params.carbon;

    %Colors
    black = [0, 0, 0];
    gray = [0.4, 0.4, 0.4];
    H2blue = [0/255, 176/255, 240/255];
    colors.H2blue = [0/255, 176/255, 240/255];
    colors.ColumbiaBlue = [155/255, 221/255, 255/255];
    colors.ChiliRed = [194/255, 24/255, 7/255];
    colors.Raspberry = [210/255, 31/255, 60/255];
    colors.MITRed = [163/255, 31/255, 52/255];
    colors.MITDarkGrey = [138/255, 139/255, 140/255];
    colors.MITLightGrey = [194/255, 192/255, 191/255];
    colors.hydrogen = [0, 176/255, 240/255];
    colors.methane = [234/255, 112/255, 13/255];
    colors.HC_linear = [192/255, 0, 0];
    colors.HC_rings = [128/255, 64/255, 64/255];
    colors.sulfur = [225/255, 173/255, 1/255];
    colors.navyblue = [0/255,0/255,128/255];
    colors.trueblue = [0/255,115/255,207/255];
    colors.elecblue = [125/255,249/255,255/255];
    colors.raspred = [179/255,68/255,108/255];
    colors.truered = [255/255,0/255,0/255];
    colors.rubyred = [155/255,17/255,30/255];

    %Specify custom color order
    


%% Evaluate mass in each level

    %Define reference state for thermodynamic calculations
    T0 = 293.15;
    p0 = 101325;
    
    %Analysis settings
    Npts_layer = 50;

    %Limit analysis to final (approx. steady) state
    ts = t;
    ts_init = t;
    t_ind = length(t);
    t = t(t_ind);
    y_t = params.Ns_m; %y(t_ind, :);
    y_t_all = y(end,:);

    %Set the initial bin to 0
    y_t(1) = 0;
    y_t_all(1) = 0;

    %Pull relevant values
    X_Ar = params.X_Ar;


    %TEMPORARY
    rep_m_ind = 20; %Representative index for calculating mean conversion and temperature

    %Pull indices
    xinds = params.xinds_m;
    zinds = params.zinds_m;
    Nz = length(unique(zinds));

    %Iterate through z-levels
    it = 1;
    d_sauter = zeros(params.Nz, 1);
    dbsout = zeros(params.Nz, params.Nms);
    mtot = zeros(1,params.Nz);
    mflux = zeros(1, params.Nz);
    Vtot = zeros(1, params.Nz); 
    alphag = zeros(1, params.Nz);
    uspf = zeros(1, params.Nz);
    mdists = zeros(params.Nz, params.Nms);
    mdists_norm = mdists;
    Tbars = zeros(params.Nz, params.Nms);
    Hbars = zeros(params.Nz, params.Nms);
    Xbars = zeros(params.Nz, params.Nms);
    nCH4s = zeros(params.Nz, params.Nms);
    nH2s = zeros(params.Nz, params.Nms);
    nCs = zeros(params.Nz, params.Nms);
    nArs = zeros(params.Nz, params.Nms);
    Tg = zeros(1, params.Nz);
    Tg_mix = zeros(1, params.Nz);
    Tdists = cell(params.Nz, 1);
    Xdists = Tdists;
    Vbs_all = Tdists;
    Hdists = zeros(params.Nz, Npts_layer);
    Hdists_smooth = Hdists;
    Tsms = zeros(params.Nz, Npts_layer);
    Xg = zeros(1, params.Nz); %Overall conversion for each z-level, calculated from the representative mass conversion - IMPROVEMENT: Calculate based on full distribution
    Xg_mix = zeros(1, params.Nz); %Mole-weighted conversion for each z-level, calculated from the total moles of CH4 and H2 in the gas phase
    %Nbs = zeros(params.Nz, params.Nms); %Storage matrix for numeric density of bubbles in each cell

    for iz = 1:params.Nz

        %Pull current z-levels
        zm = params.zms(iz);
        ub_func = params.ubs.funcs{iz};
        p = params.p_z(zm);

        %Pull the current distribution of bubble sizes
        mdists(iz,:) = y_t(zinds == iz);
        mdists_norm(iz,:) = mdists(iz,:)./max(mdists(iz,:));


        %Calculate the properties from individual bubbles
        dbs = zeros(1, params.Nms);
        mbs = zeros(1, params.Nms);
        Vbs = zeros(1, params.Nms);
        Nbs = zeros(1, params.Nms);
        Tdists_z = cell(1, params.Nms);
        Ndists_z = cell(1, params.Nms);
        
        Vdots = Vbs;
        mfluxs = mbs;
        for im = 1:params.Nms

            %Calculate current volume and diameter
            m = params.mms(im);
            ni = params.nms(im);
            X_bar = 0.5; %IMPROVEMENT - Create weighted averaging function once temp and conversion implemented
            [T_bar, Ts_m, Ns_m] = CalcTbar(y(end,:), y_t, iz, im, params);
            [X_bar, Xs_m, Ns_m_X] = CalcXbar(y(end,:), y_t, iz, im, params);
            n_Ar = ni * params.X_Ar;
            n_CH4_i = ni * (1-X_Ar); %mols - initial moles of methane in each bubble 
            n_CH4 = (1-X_bar) * n_CH4_i;
            n_H2 = (2*X_bar) * n_CH4_i;
            n_C = X_bar * n_CH4_i;
            n_gas = n_Ar + n_CH4 + n_H2;
            m_gas = n_CH4*methane.mol_mass + n_H2*hydrogen.mol_mass + n_Ar*argon.mol_mass;
            m_solid = n_C * carbon.mol_mass;
            m_b = m_gas + m_solid;
            M_bar = (m_gas)./n_gas;
            V = (n_gas * params.R * T_bar)/p;
            d = (6.*V/pi).^(1/3);

            %Calculate enthalpy for the bubble 
            H_CH4 = n_CH4 * integral(methane.Cp, T0, T_bar); %J/bubble
            H_H2 = n_H2 * integral(hydrogen.Cp, T0, T_bar);
            H_Ar = n_Ar * integral(argon.Cp, T0, T_bar);
            H_C = n_C * integral(carbon.Cp, T0, T_bar);
            H_bar = H_CH4 + H_H2 + H_Ar + H_C; %Enthalpy of the bubble based on the gas composition and temperature - I;

            %Calculate current velocity
            ub = ub_func(d);
            if isnan(ub) && d < 1E-3
                ub = ub_func(1E-3);
            end

            %Calculate number of bubbles in the ccell
            t_res_b = params.mesh.hs(iz)/ub; %s - time spent by the bubble in the cell
            ind = find(zinds == iz & xinds == im);
            N = y_t(ind) * params.Vcells(iz); %# - number of bubbles in the cell
            Vbs(im) = N * V;

            %Calculate numeric flux (simple method)
            Nflux_simple(im) = N/t_res_b;

            %Calculate a total gas flux
            if N > 0 && params.heat.active
                Nfluxs(im) = 0;
                Ns_m_frac = (Ns_m)./sum(Ns_m);
                for it = 1:length(Ns_m)
                    N_m_frac = Ns_m_frac(it);
                    if N_m_frac > 0

                        %Pull temperature and conversion
                        T = Ts_m(it);
                        X = Xs_m(it);

                        %Calculate moles of gas
                        n_Ar = ni * params.X_Ar;
                        n_CH4_i = ni * (1-X_Ar); %mols - initial moles of methane in each bubble 
                        n_CH4 = (1-X) * n_CH4_i;
                        n_H2 = (2*X) * n_CH4_i;
                        n_C = X * n_CH4_i;
                        n_gas_T = n_Ar + n_CH4 + n_H2;

                        %Calculate volume and diameter
                        V_T = (n_gas_T * params.R * T)/p;

                        %Calculate velocity and residence time
                        d_T = (6.*V_T/pi).^(1/3);
                        ub_T = ub_func(d_T);
                        t_res_b_T = params.mesh.hs(iz)/ub_T %s - residence time

                        %Add to the numerical flux
                        Nfluxs(im) = Nfluxs(im) + N*N_m_frac/t_res_b_T;


                    end
                    
                end

                %Nfluxs(im) = 1;
                %mfluxs(im) = 1;

            else
                Nfluxs(im) = Nflux_simple(im); %Number of bubbles entering/exiting per second
                
            end


            
            

            %Pull numeric densities for this representative mass
            subinds = find(zinds == iz & xinds == im);
            subNs = y_t(subinds); %All subbins for this mass



            %Calculate mass in cell and mass flux in/out
            dbs(im) = d;
            mbs(im) = N*m_b;
            Nbs(im) = subNs;
            mfluxs(im) = Nfluxs(im) * m_b; %kg/s
            Vdots(im) = Nfluxs(im) * V;
            Tbars(iz,im) = T_bar;
            Hbars(iz,im) = H_bar;
            Xbars(iz,im) = X_bar;
            nCH4s(iz,im) = n_CH4;
            nH2s(iz,im) = n_H2;
            nCs(iz,im) = n_C;
            nArs(iz,im) = n_Ar;
            Tdists_z{im} = Ts_m;
            Ndists_z{im} = Ns_m;
            Xdists_z{im} = Xs_m;



            it = it + 1;

        end

        %Log results
        Tdists{iz} = Tdists_z;
        Ndists{iz} = Ndists_z;
        Vbs_all{iz} = Vbs;

        %Calculate sauter mean diameter
        d_sauter(iz) = sum(Nbs .* dbs.^3)./sum(Nbs .* dbs.^2);


        %Calculate fraction of each size i
        dbsout(iz,:) = dbs;
        Vdot(iz) = sum(Vdots);
        Vtot(iz) = sum(Vbs);
        Tg(iz) = Tbars(iz,:) * [Nbs./sum(Nbs)]'; %Weighted average gas temperature for each z-level - UPDATE TO CALCULATE MIXING CUP TEMP
        Xg(iz) = Xbars(iz,:) * [Nbs./sum(Nbs)]';
            


        %Calculate mean composition of 
        nCH4g(iz) = nCH4s(iz,:) * Nbs'; %mol/m3 - Total number of moles of CH4 in the gas phase for each z-level, calculated from the representative mass composition and numeric density of bubbles in each cell
        nH2g(iz) = nH2s(iz,:) * Nbs';
        nCg(iz) = nCs(iz,:) * Nbs';
        nArg(iz) = nArs(iz,:) * Nbs';

        %Calculate mole weighted conversion
        Xg_mix(iz) = (nH2g(iz)/2)./(nCH4g(iz) + nH2g(iz)/2); %Mole-weighted conversion for each z-level, calculated from the total moles of CH4 and H2 in the gas phase

        %Calculate mixing cup temperature of the gas for this later
        Hg(iz) = Hbars(iz,:) * Nbs'; %J/m3 - total enthalpy of the layer 
        Tg_func = @(T) nCH4g(iz) * integral(methane.Cp, T0, T) + nH2g(iz) * integral(hydrogen.Cp, T0, T)...
                     + nCg(iz) * integral(carbon.Cp, T0, T) + nArg(iz) * integral(argon.Cp, T0, T) - Hg(iz);
        Tg_mix(iz) = fzero(Tg_func, Tg(iz)); %K - mixing cup temperature of the gas, calculated from the total enthalpy and composition of the gas in the layer

        %Calculate total temperature distribution for this layer
        
        if params.heat.active
            
            T_min_layer = min(Tdists_z{end}); %K - minimum temperature of all bubbles in this layer (found in the largest bubble)
            T_max_layer = params.T_liq; %K - maximum temperature of all bubbles in this layer (likely liquid temperature)
            Tsb_layer = linspace(T_min_layer, T_max_layer, Npts_layer+1); %K - a normal temperature distribution for calculating the total temperature distribution
            Tsm_layer = Tsb_layer(1:end-1) + diff(Tsb_layer)./2;
            Hdist_layer = zeros(1, Npts_layer); %Storage for the total enthalpy distribution for this layer, calculated by allocating the enthalpy of each bubble into the appropriate temperature bins based on their individual temperature distributions
            for im = 1:params.Nms
                Vm = params.Vms(im); %m3 - volume for weighting
                Tdist = Tdists_z{im};
                Ndist_T = Ndists_z{im};
                Xdist = Xdists_z{im};
                Hdist = zeros(1, length(Tdist));

                %Calculate enthalpy distribution
                for it = 1:length(Tdist)

                    %Pull temperature distribution
                    T = Tdist(it);
                    

                    %Determine mean composition of the gas at this moment -
                    %TO BE IMPLEMENTED
                    X_Ar = 0;
                    X_CH4 = 0;
                    X_H2 = 2/3; %Mole fraction of hydrogen
                    X_C = 1/3; %Solid carbon enthalpy needs to be accounted for

                    %Calculate number of moles in the control volume 
                    nf_ratio = 1;

                    %Calculate enthalpy at this temperature                 
                    Hdist(it) = X_Ar * integral(argon.Cp,params.T0,T)...
                        + X_CH4 * integral(methane.Cp,params.T0, T)...
                        + X_C * integral(carbon.Cp,params.T0, T)...
                        + X_H2 * integral(hydrogen.Cp,params.T0, T);

                end

                %Calculated weighted enthalpy distribution
                Hdist = Hdist .* Ndist_T .* Vm;

                %Allocate to bins
                for itl = 1:Npts_layer

                    %Define bounds
                    Tu = Tsb_layer(itl+1);
                    Tl = Tsb_layer(itl);

                    %Determine if any bubbles from this mass fall within the bounds
                    inds_within = find(Tdist >= Tl & Tdist < Tu);

                    %If so, add their enthalpy to the total enthalpy distribution for this layer
                    if ~isempty(inds_within)
                        Hdist_layer(itl) = Hdist_layer(itl) + sum(Hdist(inds_within));
                    end


                end

                x = 1;

            
                
            end

            %Smooth distribution - CDF to PDF
            Hdist_cdf = cumsum(Hdist_layer);
            Hdist_cdf_spline = spline(Tsm_layer, Hdist_cdf);
            Hdist_spline_coeffs = Hdist_cdf_spline.coefs .* repmat([3,2,1,0], Hdist_cdf_spline.pieces, 1);
            Hdist_spline_coeffs = [zeros(Hdist_cdf_spline.pieces, 1), Hdist_spline_coeffs(:,1:3)];
            Hdist_spline = mkpp(Hdist_cdf_spline.breaks, Hdist_spline_coeffs);
            Hdist_spline_vals = ppval(Hdist_spline, Tsm_layer);
            Hdist_ks = ksdensity(Hdist_layer);

            %Smooth distribution - lowess smoothing
            span = round(0.25*Npts_layer);
            Hdist_smooth = smooth(Tsm_layer, Hdist_layer, span, 'lowess');%sgolayfilt(Hdist_layer,3,15); %fit(Tsm_layer', Hdist_layer', 'smoothingspline')
            
            
            %Hdist_movmean = movmean(Hdist_layer, 21);


            %Log temperature distributions
            Tsms(iz,:) = Tsm_layer;
            Hdists(iz, :) = Hdist_layer./sum(Hdist_layer);
            Hdists_smooth(iz,:) = Hdist_smooth./sum(Hdist_smooth);

        end

        


        %Calculate mean diameter
        dbar(iz) = sum(Nbs .* dbs.^3)./sum(Nbs .* dbs.^2); %m - Sauter mean diameter for each z-level, calculated from the numeric density of bubbles and their diameters

        %Calculate gas holdup, superficial velocity, and average bubble velocity for this layer
        alphag(iz) = Vtot(iz)/params.Vcells(iz);
        uspf(iz) = Vdot(iz)/params.reactor.Ac; %m/s - superficial velocity for each z-level, calculated from volumetric flow rate and cross-sectional area
        ub_bar(iz) = uspf(iz)/alphag(iz); %m/s - average bubble velocity for each z-level, calculated from superficial velocity and holdup
        mtot(iz) = sum(mbs);
        mflux(iz) = sum(mfluxs);

        x = 1;
    end

    %Calculate mass relative change
    mrel = mflux./mflux(1);
    
    %Log results
    results.t = ts_init; results.y = y;
    results.dsauter = d_sauter;
    results.dbsout = dbsout;
    results.Tg_mix = Tg_mix;
    results.Xg_mix = Xg_mix;
    results.Vdot = Vdot;
    results.Vtot = Vtot;
    results.alphag = alphag;
    results.alphag_avg = sum(Vtot)/(params.reactor.H*pi*params.reactor.R^2);
    results.uspf = uspf;
    results.mtot = mtot;
    results.mflux = mflux;
    results.mdists = mdists;
    results.mdists_norm = mdists_norm;
    results.params = params;
    




%% Normalize based on integral
    
    mdists_normint = zeros(size(mdists_norm));

    for iz = 1:params.Nz
        mdist = mdists(iz,:);
        dbs = dbsout(iz,:) .* 1000; %Convert to mm
        mdist_int = trapz(dbs, mdist);
        mdists_normint(iz,:) = mdist./mdist_int;

    end

%% 


%% Plots

if params.sol.postprocess

    %Handle initial zero mass 
    if isnan(results.mflux)
        results.mflux(1) = results.mflux(2);
    end

    %Plot mass flow rate continuity
    figure();

    subplot(1,2,1)
    plot(1:params.Nz, 100.*results.mflux/results.mflux(1), 'k-', 'LineWidth', lw);
    xlabel('Z-cell'); ylabel('Mass flux');
    grid on; grid minor; axis square;
    title('Normalized Mass Flux');
    set(gca, 'FontSize', fs);

    %Plot for total mass timeseries
    %figure();
    subplot(1,2,2);
    m_inds_end = min([length(params.m_total), length(ts)]);
    m_inds = 1:m_inds_end; %1:length(ts(ts>=0.1));
    plot(ts(m_inds), 100.*params.m_total(m_inds)./(params.m_total(m_inds(end))), 'k-', 'LineWidth', lw);
    xlabel('Sim. Time (s)'); ylabel('Normalized Mass');
    title('Total Mass in System')
    grid on; grid minor; axis square;
    set(gca, 'FontSize', fs)

    figure()
    subplot(1,2,1);
    plot(1:params.Nz, 100.*results.uspf, 'LineWidth', lw);
    xlabel('Z-cell'); ylabel('Superficial Vel (cm/s)');
    grid on; grid minor; axis square;
    title('Gas Superficial Velocity');
    set(gca, 'FontSize', fs);

    %Plot for total mass timeseries
    % %figure();
    % m_inds_end = min([length(params.m_total), length(ts)]);
    % m_inds = 1:m_inds_end; %1:length(ts(ts>=0.1));
    % plot(ts(m_inds), 100.*params.m_total(m_inds)./(params.m_total(m_inds(end))), 'k-', 'LineWidth', lw);
    % xlabel('Sim. Time (s)'); ylabel('Normalized Mass');
    % title('Total Mass in System')
    % grid on; grid minor; axis square;
    % set(gca, 'FontSize', fs)



    %Plot normalized bubble size distributions
    figure();
    colororder('reef');

    subplot(1,2,1)
    semilogx(params.mms, mdists', 'LineWidth', lw)
    grid on; grid minor; axis square;
    title('Numeric density ')

    subplot(1,2,2)
    semilogx(params.mms, mdists_norm', 'LineWidth', lw)
    grid on; grid minor; axis square;
    title('Normalized Bubble Size Distributions')

    %Surface plot
    if ~params.sol.single_layer

        %Surface plot - normalized
        figure()
        mesh = params.mesh;
        zsc = mesh.volcell_cents(:,2);
        zsc = repmat(zsc, 1, params.Nms);
        xsc = repmat(params.mmesh.xsc, mesh.N_cells,1);
        Yf = mdists_norm;
        %Yf = reshape(Yf,params.Nms, mesh.N_cells); Yf = Yf';
        surf(100.*zsc, 100.*xsc, Yf);
        colormap jet
        ylabel('Bubble Mass (\mug)')
        xlabel('Vertical Position (cm)')
        zlabel('Normalized BSD');
        set(gca, 'YScale', 'log');
        axis square;
    

        %Surface plot - unnormalized
        figure()
        mesh = params.mesh;
        zsc = mesh.volcell_cents(:,2);
        zsc = repmat(zsc, 1, params.Nms);
        xsc = repmat(params.mmesh.xsc, mesh.N_cells,1);
        Yf = mdists;
        %Yf = reshape(Yf,params.Nms, mesh.N_cells); Yf = Yf';
        surf(100.*zsc, 100.*xsc, Yf);
        colormap jet
        ylabel('Bubble Mass (\mug)')
        xlabel('Vertical Position (cm)')
        zlabel('Normalized BSD');
        set(gca, 'YScale', 'log');
        axis square

    end

    %Diameter plot
    figure(); %figure('WindowStyle', 'docked');
    colororder('reef');

    subplot(1,2,1)
    plot(1000.*dbsout(end,:), mdists', 'LineWidth', lw)
    grid on; grid minor; axis square;
    xlabel('Diameter (m)');
    title('Numeric density ')

    subplot(1,2,2)
    plot(1000.*dbsout(end,:), mdists_normint', 'LineWidth', lw)
    grid on; grid minor; axis square;
    xlabel('Diameter (m)')
    title('Normalized Bubble Size Distributions')

    %Single layer timeseries plot
    dt = mode(diff(ts));
    N_ts = length(ts);
    N_plots = 10;
    dN = round((N_ts-1)/(N_plots-1));
    
    if params.sol.single_layer

        figure();



        %Define colors
        newcolors = [colors.navyblue;
                     colors.trueblue;
                     colors.H2blue;
                     colors.ColumbiaBlue;
                     colors.MITDarkGrey;
                     colors.truered;
                     colors.MITRed;
                     colors.ChiliRed];
        colororder(newcolors)

        
        in = 1;
        ts_plot = [];
        tstrs = {};
        for i = 1:N_plots
            ts_plot(end+1) = ts(in);
            tstrs{i} = sprintf('t = %0.4f s', ts(in));
            y_plot = y(in,:);
            subplot(1,2,1);
            plot(dbsout(end,:), y_plot, 'LineWidth', lw); hold on;
            subplot(1,2,2);
            semilogy(dbsout(end,:), y_plot, 'LineWidth', lw); hold on;
            

    
            in = in + dN;
            if in > length(ts)
                in = length(ts);
            end

        end

        subplot(1,2,1);
        xlabel('Diameter (m)'); ylabel('Numeric Density')
        grid on; grid minor; axis square;
        legend(tstrs{1}, tstrs{2}, tstrs{3}, tstrs{4}, tstrs{5}, tstrs{6}, tstrs{7}, tstrs{8});
        set(gca, 'FontSize', fs, 'FontWeight', 'bold');
        title('Single Layer - Timeseries')
        ylim([0, 18E4]);

        subplot(1,2,2);
        xlabel('Diameter (m)'); ylabel('Numeric Density')
        grid on; grid minor; axis square;
        legend(tstrs{1}, tstrs{2}, tstrs{3}, tstrs{4}, tstrs{5}, tstrs{6}, tstrs{7}, tstrs{8});
        set(gca, 'FontSize', fs, 'FontWeight', 'bold');
        ylim([1E-3, 10000])
        title('Single Layer - Timeseries')
        
    end

    %Probability density (in terms of volume) plot
    figure();
    subplot(1,2,1);
    V_total = sum(Vbs);
    V_total_i = sum(Vbs_all{1});
    plot(1000*dbsout(1,:), Vbs_all{1}./V_total_i, 'LineWidth', lw, 'Color', colors.MITDarkGrey); hold on;
    plot(dbs, Vbs./V_total, 'k-', 'LineWidth', lw, 'color', colors.H2blue);
    grid on; grid minor; axis square;

    xlabel('Diameter (m)'); ylabel('PDF (volume) of d_b');
    legend('Initial', 'Final', 'location', 'southeast');
    title('PDF (volume/mass) of d_b')
    set(gca, 'FontSize', fs);

%% Plot mean quantities versus height

    if ~params.sol.single_layer

        %Gas holdup and superficial velocity versus height
        figure('Name', 'Gas Holdup and Superficial Velocity vs Height');
    
        subplot(1,2,1)
        plot(100.*zsc(:,1), alphag, 'LineWidth', lw, 'Color', colors.H2blue); hold on;
        xlabel('Vertical Position (cm)'); ylabel('Gas Holdup');
        grid on; grid minor; axis square;
        title('Gas Holdup vs Height');
        set(gca, 'FontSize', fs);
    
        subplot(1,2,2);
        plot(100.*zsc(:,1), uspf, 'LineWidth', lw, 'Color', colors.ChiliRed); hold on;
        xlabel('Vertical Position (cm)'); ylabel('Superficial Velocity (cm/s)');
        grid on; grid minor; axis square;
        title('Superficial Velocity vs Height');
        set(gca, 'FontSize', fs);

        %Plot mean bubble velocity and diameter 
        figure('Name', 'Mean Bubble Velocity and Diameter vs Height');
    
        subplot(1,2,1);
        plot(100.*zsc(:,1), ub_bar, 'LineWidth', lw, 'Color', colors.MITRed); hold on;
        xlabel('Vertical Position (cm)'); ylabel('Average Bubble Velocity (cm/s)');
        grid on; grid minor; axis square;
        title('Average Bubble Velocity vs Height');
        set(gca, 'FontSize', fs);

        subplot(1,2,2);
        plot(100.*zsc(:,1), dbar.*1000, 'LineWidth', lw, 'Color', colors.MITDarkGrey); hold on;
        xlabel('Vertical Position (cm)'); ylabel('Sauter Mean Diameter (mm)');
        grid on; grid minor; axis square;
        title('Sauter Mean Diameter vs Height');
        set(gca, 'FontSize', fs);

        %Plot gas composition versus height
        figure()

        subplot(1,2,1);
        plot(100.*zsc(:,1), nCH4g, 'LineWidth', lw, 'Color', colors.methane); hold on;
        plot(100.*zsc(:,1), nH2g, 'LineWidth', lw, 'Color', colors.hydrogen); hold on;
        plot(100.*zsc(:,1), nCg, 'LineWidth', lw, 'Color', colors.raspred); hold on;
        plot(100.*zsc(:,1), nArg, 'LineWidth', lw, 'Color', colors.navyblue); hold on;
        xlabel('Vertical Position (cm)'); ylabel('Moles of Species (mol/m3)');
        grid on; grid minor; axis square;
        title('Gas Composition vs Height');
        legend('CH_4', 'H_2', 'C', 'Ar');

        subplot(1,2,2)
        %plot()

    end



%% Plot temperature and conversion versus height

    %Create figure
    figure('Name', 'Temperature and Conversion vs Height');

    %Plot settings
    dim = 3;

    %Specify masses to sample
    msamps = [1, 5:5:params.Nms];

    %Specify colors
    
    if params.heat.active
    
        %For a sample of representative masses, plot Tbar and Xbar versus height
        for i = 1:length(msamps)

            %Define mass index
            im = msamps(i);

            %Plot mean temperature
            subplot(1,2,1)
            plot(100.*zsc(:,1), Tbars(:,im)-273.15, 'LineWidth', lw); hold on;
            %plot(100.*zsc(:,1), Tg_mix(:), 'LineWidth', lw); hold on;
            %plot(100.*zsc(:,1), Tg(:), 'LineWidth', lw, 'Color', colors.MITDarkGrey); hold on;
            xlabel('Vertical Position (cm)'); ylabel('Tbar (C)');
            grid on; grid minor; axis square;
            title('Tbar vs Height');
            set(gca, 'FontSize', fs);
    
            %Plot mean conversion
            subplot(1,2,2)
            plot(100.*zsc(:,1), 100.*Xbars(:,im), 'LineWidth', lw); hold on;
            xlabel('Vertical Position (cm)'); ylabel('Xbar');
            grid on; grid minor; axis square;
            title('Xbar vs Height');
            set(gca, 'FontSize', fs);

            %Specify legend entires
            legend_entires{i} = sprintf('im = %d; d = %0.4f m', im, dbsout(end,im));


        end

        %Add legend entires to both plots
        subplot(1,2,1)
        legend(legend_entires, 'location', 'southeast');
        subplot(1,2,2)
        legend(legend_entires, 'location', 'southeast');

        %Plot the weighted average temperature (Tg) and conversion (Xg) versus height
        figure('Name', 'Temp. and Conv. vs Height');
        subplot(1,2,1)
        plot(100.*zsc(:,1), Tg_mix-273.15, 'r-', 'LineWidth', lw, 'Color', colors.ChiliRed); hold on;
        plot(100.*zsc(:,1), Tg-273.15, 'r-', 'LineWidth', lw, 'Color', colors.hydrogen);
        xlabel('Vertical Position (cm)'); ylabel('Tg (C)');
        grid on; grid minor; axis square;
        title('Gas Temp. vs Height');
        legend('Tg - Mixing Cup Temp', 'Tg - Weighted Avg Temp', 'location', 'southeast');
        set(gca, 'FontSize', fs);

        subplot(1,2,2);
        plot(100.*zsc(:,1), 100.*Xg, 'k-', 'LineWidth', lw, 'Color', colors.hydrogen); hold on;
        plot(100.*zsc(:,1), 100.*Xg_mix, 'r-', 'LineWidth', lw, 'Color', colors.ChiliRed); hold on;
        xlabel('Vertical Position (cm)'); ylabel('Xg');
        grid on; grid minor; axis square;
        title('Gas Conversion vs Height');
        legend('Xg - Weighted Avg', 'Xg - Mole Weighted', 'location', 'southeast');
        set(gca, 'FontSize', fs);


    end


%% Plot temperature distribution for each height

    if params.heat.active

        figure('Name', 'Temperature Distributions vs Height');

        %Plot the final distribution
        subplot(1,2,1);
        bar(Tsms(end,:)-273.15, Hdists(end,:), 'FaceColor', colors.ChiliRed, 'EdgeColor', colors.ChiliRed);
        xlabel('Temperature (C)'); ylabel('PDF');
        grid on; grid minor; axis square;
        title(sprintf('Final Temp. Distribution (Z = %0.2f cm)', 100.*zsc(end,1)));

        %Define 10 colors in single scale
        Nz = params.Nz;
        colors_heat = [linspace(0,1,Nz)', 0.2*ones(Nz,1), flipud(linspace(0,1,Nz)')];

        %Plot each layer 
        subplot(1,2,2);
        for iz = 1:params.Nz
        
            plot(Tsms(iz,:)-273.15, Hdists_smooth(iz,:), 'k-', 'Marker', '.', 'MarkerSize', 8, 'Color', colors_heat(iz,:)); hold on;
            xlabel('Temperature (C)'); ylabel('PDF');
            grid on; grid minor; axis square;
            title(sprintf('Z = %0.2f cm', 100.*zsc(iz,1)));

            %Create legend entry
            legend_entries{iz} = sprintf('Z = %0.2f cm', 100.*zsc(iz,1));

        end

        %Add legend
        legend(legend_entries, 'location', 'northeast');
        


    end


%% Summary Plot - Multilayer 

if ~params.sol.single_layer
    %Create figure
    figure('Name', 'Summary Plot - Multilayer', 'units', 'normalized', 'outerposition', [0, 0, 10/16, 1]);

    %Plot 3D graph - numeric density versus diameter and height
    subplot(2,2,1)
    mesh = params.mesh;
    zsc = mesh.volcell_cents(:,2);
    zsc = repmat(zsc, 1, params.Nms);
    xsc = repmat(dbsout(end,:), mesh.N_cells,1);
    Yf = mdists;
    %Yf = reshape(Yf,params.Nms, mesh.N_cells); Yf = Yf';
    surf(100.*zsc, 100.*xsc, Yf);
    colormap jet
    ylabel('Bubble Diameter (cm)')
    xlabel('Vertical Position (cm)')
    zlabel('Normalized BSD');
    title('BSD vs Height')
    set(gca, 'YScale', 'log');
    axis square


    %Plot final bubble size distribution numeric density (yyaxis left) and volume pdf (yyaxis right) for the final z-level 
    subplot(2,2,2)
    %yyaxis left 
    %mdist_in = params.N_dot_o;
    plot(dbsout(1,:), mdists(1,:), 'LineWidth', lw, 'Color', colors.MITDarkGrey); hold on;
    plot(dbsout(end,:), mdists(end,:), 'LineWidth', lw, 'Color', colors.hydrogen);
    ylabel('Numeric Density (1/m^3)'); xlabel('Bubble Diameter (m)');
    title('Bubble Size Distribution (BSD)');
    %ylim([1, 10*max(mdists(1,:))]);
    grid on; grid minor; axis square;
    legend('Initial', 'Final', 'location', 'southeast');

    axes('Position', [0.71, 0.75, 0.15, 0.15]);
    box on;
    plot(dbsout(1,:), mdists(1,:), 'LineWidth', lw, 'Color', colors.MITDarkGrey); hold on;
    plot(dbsout(end,:), mdists(end,:), 'LineWidth', lw, 'Color', colors.hydrogen);
    ylim([0, 0.1*max(mdists(end,:))])

    %Decide what to plot in the 3 and 4 sections 
    if params.heat.active %Include temp/conv

        %Plot mixing cup temperature (yyaxis left) and conversion (yyaxis right) versus height 
        subplot(2,2,3)
        yyaxis left
        plot(100.*zsc(2:end,1), Tg_mix(2:end)-273.15, 'r-', 'LineWidth', lw, 'Color', colors.ChiliRed); hold on;
        plot(100.*zsc(2:end,1), Tg(2:end)-273.15, 'r-.', 'LineWidth', lw, 'Color', colors.ChiliRed); hold on;   
        ylabel('Temperature (C)')
        grid on; grid minor; axis square;
        set(gca,'YColor', colors.ChiliRed);
        yyaxis right
        plot(100.*zsc(2:end,1), 100.*Xg_mix(2:end), 'r-', 'LineWidth', lw, 'Color', colors.hydrogen); hold on;
        plot(100.*zsc(2:end,1), 100.*Xg(2:end), 'k-.', 'LineWidth', lw, 'Color', colors.hydrogen); hold on;
        set(gca,'YColor', colors.hydrogen);
        ylabel('Conversion (%)')
        xlabel('Vertical Position (cm)')
        title('Temp. and Conv. vs Height');
        legend('Gas Temp. - Mixing Cup', 'Gas Temp. - Mean', 'Gas Conv. - Mole Weighted', 'Gas Conv. - Mean', 'location', 'southeast', 'BackgroundAlpha', 0.75);


        %Plot final temperature distribution for the final z-level
        subplot(2,2,4);
        bar(Tsms(end,:)-273.15, Hdists(end,:), 'FaceColor', colors.ChiliRed, 'EdgeColor', colors.ChiliRed, 'FaceAlpha', 0.5); hold on;
        plot(Tsms(end,:)-273.15, Hdists_smooth(end,:), 'k-.', 'LineWidth', lw, 'Color', colors.ChiliRed); hold on;
        xlabel('Temperature (C)'); ylabel('PDF');
        grid on; grid minor; axis square;
        title(sprintf('Final Temp. Distribution'));
        legend('Temp. Dist. - Raw', 'Temp. Dist. - Smoothed', 'location', 'northwest');
        ylim([-0.01, max(Hdists(end,:) + 0.03)])

    else %Don't include temp/conv
        x  =1;


    end


end

%% Summary plot - Single-layer

if params.sol.single_layer

    figure('Name', 'Summary Plot - Single Layer');

    %Plot the initial and final bubble size distributions
    subplot(2,2,1);
    plot(dbsout(end,:), mdists(1,:), 'LineWidth', lw, 'Color', colors.H2blue); hold on;
    plot(dbsout(end,:), mdists(end,:), 'LineWidth', lw, 'Color', colors.ChiliRed); hold on;



end


end

end