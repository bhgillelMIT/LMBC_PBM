function results = PBM_postprocess(t, y, params)
% A function to process the outputs from the PBM. Analyzes the bubble size
% distribution, volumetric and mass flow rates, and other physical aspects
% of the simulation.


%% Setup

    %Handle no inputs
    if nargin < 1
        output = load('Data/Solutions/PBM_output_16-Dec-2025_15-58-33.mat');
        output = output.output;
        t = output.T;
        y = output.Y;
        params = output.params;
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

    %Limit analysis to final (approx. steady) state
    ts = t;
    t_ind = length(t);
    t = t(t_ind);
    y_t = y(t_ind, :);

    %Pull relevant values
    X_Ar = params.X_Ar;


    %TEMPORARY
    rep_m_ind = 20; %Representative index for calculating mean conversion and temperature

    %Pull indices
    xinds = params.xinds;
    zinds = params.zinds;
    Nz = length(unique(zinds));

    %Iterate through z-levels
    it = 1;
    dbsout = zeros(params.Nz, params.Nms);
    mtot = zeros(1,params.Nz);
    mflux = zeros(1, params.Nz);
    Vtot = zeros(1, params.Nz); 
    alphag = zeros(1, params.Nz);
    uspf = zeros(1, params.Nz);
    mdists = zeros(params.Nz, params.Nms);
    mdists_norm = mdists;

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
        Vdots = Vbs;
        mfluxs = mbs;
        for im = 1:params.Nms

            %Calculate current volume and diameter
            m = params.mms(im);
            ni = params.nms(im);
            X_bar = 0.5; %IMPROVEMENT - Create weighted averaging function once temp and conversion implemented
            T_bar = CalcTbar(y_t, iz, im, params);
            X_bar = CalcXbar(y_t, iz, im, params);
            n_Ar = ni * params.X_Ar;
            n_CH4_i = ni * (1-X_Ar);
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

            %Calculate current velocity
            ub = ub_func(d);

            %Calculate number of bubbles in the ccell
            t_res_b = params.mesh.hs(iz)/ub; %s - time spent by the bubble in the cell
            ind = find(zinds == iz & xinds == im);
            N = y_t(ind) * params.Vcells(iz); %# - number of bubbles in the cell
            Vbs(im) = N * V;



            %Calculate mass in cell and mass flux in/out
            dbs(im) = d;
            mbs(im) = N*m_b;
            Nfluxs(im) = N/t_res_b; %Number of bubbles entering/exiting per second
            mfluxs(im) = Nfluxs(im) * m_b; %kg/s
            Vdots(im) = Nfluxs(im) * V;

            it = it + 1;

        end
        dbsout(iz,:) = dbs;
        Vdot(iz) = sum(Vdots);
        Vtot(iz) = sum(Vbs);
        alphag(iz) = Vtot(iz)/params.Vcells(iz);
        uspf(iz) = Vdot(iz)/params.reactor.Ac;
        mtot(iz) = sum(mbs);
        mflux(iz) = sum(mfluxs);

        x = 1;
    end

    %Calculate mass relative change
    mrel = mflux./mflux(1);
    


    %Log results
    results.t = t; results.y = y;
    results.dbsout = dbsout;
    results.Vdot = Vdot;
    results.Vtot = Vtot;
    results.alphag = alphag;
    results.uspf = uspf;
    results.mtot = mtot;
    results.mflux = mflux;
    results.mdists = mdists;
    results.mdists_norm = mdists_norm;
    


%% Normalize based on integral
    
    mdists_normint = zeros(size(mdists_norm));

    for iz = 1:params.Nz
        mdist = mdists(iz,:);
        dbs = dbsout(iz,:) .* 1000; %Convert to mm
        mdist_int = trapz(dbs, mdist);
        mdists_normint(iz,:) = mdist./mdist_int;

    end

%% Plots

    %Plot mass flow rate continuity
    figure();

    subplot(1,2,1)
    plot(1:params.Nz, 100.*results.mflux/results.mflux(1), 'LineWidth', lw);
    xlabel('Z-cell'); ylabel('Mass flux');
    grid on; grid minor; axis square;
    title('Normalized Mass Flux');
    set(gca, 'FontSize', fs);

    subplot(1,2,2);
    plot(1:params.Nz, 100.*results.uspf, 'LineWidth', lw);
    xlabel('Z-cell'); ylabel('Superficial Vel (cm/s)');
    grid on; grid minor; axis square;
    title('Gas Superficial Velocity');
    set(gca, 'FontSize', fs);

    %Plot for total mass timeseries
    figure();
    m_inds = 2:1:length(find(params.m_total > 0));
    plot(ts, 100.*params.m_total(m_inds)./(params.m_total(1)), 'k-', 'LineWidth', lw);
    xlabel('Sim. Time (s)'); ylabel('Normalized Mass');
    title('Total Mass in System')
    grid on; grid minor; axis square;



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

    end

    %Diameter plot
    figure();
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
    N_plots = 6;
    dN = round(N_ts/N_plots);
    
    if params.sol.single_layer
        figure();
        in = 1;
        ts_plot = [];
        for i = 1:N_plots
            ts_plot(end+1) = ts(in);
            y_plot = y(in,:);
            plot(dbsout(end,:), y_plot, 'LineWidth', lw); hold on;
            

    
            in = in + dN;

        end
        xlabel('Diameter (m)'); ylabel('Numeric Density')
        

    end


end