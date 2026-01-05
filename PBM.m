function PBM()

%% Setup

    %Clean up
    close all

    %Debug
    debug = true;

    %Plot settings
    lw = 2;
    fs = 18;

    %Material Properties
    addpath('Material Properties/');
    load tin.mat
    load hydrogen.mat
    load carbon.mat
    load methane.mat

%% Inputs

    %Animation settings
    animate = true;
    saveanimation = true;
    fname = 'PBM_gif_v2.gif';
    gifdelay = 0.05;
    time_dilation = 1;

    %Simulation
    t_end = 3.5;
    rel_tol = 1E-6;
    abs_tol = 1E-6;
    mesh_type = 'Hybrid'; %Type of mesh to use. Options: 'Linear', 'Geometric', 'Hybrid'
    mesh_hybrid_cells = 10;
    mesh_hybrid_frac = 0.2;


    %Geometry
    reactor.D = 0.05;
    reactor.R = reactor.D./2;
    reactor.H = 1;
    reactor.Ac = pi .* reactor.R^2;

    %Physical constants
    R = 8.3145;
    g = 9.81;

    %Operating parameters
    reactor.T = 1300 + 273.15;
    reactor.T_orifice = reactor.T;
    reactor.T_min = 250 + 273.15;
    reactor.H = 1;
    reactor.X_Ar = 0.9;
    reactor.u_spf_orifice = 0.05; %m/s
    p_surf = 101325;
    rho_L = 6500;
    M_gas = 0.016;

    %Calculate reactor floew rate
    reactor.V_dot = reactor.u_spf_orifice .* reactor.Ac; %m3/s
    reactor.V_dot_Lmin = 1000 .* reactor.V_dot .* 60;


    %Discretization
    Nr = 1;
    Nz = 30; 
    Nms = 30;
    Nrs = Nms; % # of discrete volumes considered
    
    r_min = 1 .* 0.001; %m - smallest size of bubble to consider
    r_max = 15 .* 0.001;
    dt = 0.01; %initial step size, solver is adaptive


    %Initial condition
    u_i = 0;
    T_i = 573.15;


    %Boundary conditions
    mu_i = 0.0075;
    std_i = 0.001;

    %Calculate gas density at the orifice
    p_orifice = p_surf + rho_L * g * reactor.H;
    rho_gas_orifice = (p_orifice * M_gas)./(R * reactor.T_orifice);

    %Calculate mass - case with mass as size dimension
    V_min = (4/3) .* pi .* r_min^3;
    V_max = (4/3) .* pi .* r_max^3;
    m_min = V_min * rho_gas_orifice;
    m_max = V_max * rho_gas_orifice;
    switch mesh_type
        case 'Geometric'
            m_rat = (m_max/m_min).^(1./Nms);
            mbs = m_min .* m_rat.^(0:1:Nms);
        case 'Linear'
            mbs = linspace(m_min, m_max, Nms+1);
        case 'Hybrid'
            m_mid = m_min + (m_max - m_min) .* mesh_hybrid_frac;
            mbs_small = linspace(m_min, m_mid, mesh_hybrid_cells);
            m_rat = (m_max/mbs_small(end)).^(1./(Nms-mesh_hybrid_cells));
            mbs_large = mbs_small(end) .* m_rat.^(1:(Nms-mesh_hybrid_cells+1));
            mbs = [mbs_small, mbs_large];

    end
    mms = (mbs(1:end-1) + mbs(2:end))./2; %(mbs(1:end-1) + diff(mbs))./2;
    rbs = ((3*mbs*R*reactor.T_orifice)./(4*pi*M_gas*p_orifice)).^(1/3);
    rms = ((3*mms*R*reactor.T_orifice)./(4*pi*M_gas*p_orifice)).^(1/3);
    Vms = (mms .* R .* reactor.T_orifice)./(p_orifice .* M_gas);

    %Calculate radii - case with radius as size dimension
    rb_brackets = linspace(r_min, r_max, Nrs+1); %logspace(log10(r_min), log10(r_max), Nrs+1);
    rb_mids = rb_brackets(1:end-1) + diff(rb_brackets)./2;

    %Define initial distribution
    V_cell_o = reactor.Ac * (reactor.H/Nz); %m^3 - Volume of boundary cell(s)
    Fb_i = (1./(std_i .* sqrt(2.*pi))) .* exp(-0.5 .* ((rms - mu_i).^2)./(std_i.^2)); % # of bubbles/m3
    Vb_i = Vms .* Fb_i .* V_cell_o;
    Vb_i_sum = sum(Vb_i);
    N_bubbles_s = reactor.V_dot/Vb_i_sum; %Number of bubbles being generated - multiply by distribution
    Fb_i = N_bubbles_s * Fb_i;
    

    %Define PBM data storage array
    PBM_val.eg = 0; % Gas fraction
    PBM_val.uz_L = 0; % Liquid vertical velocity
    PBM_val.Fbs = zeros(1,Nrs); % Bubble volume - Discrete quantities
    PBM_val.Ts.mu = T_i; % Temperature - Mean
    PBM_val.Ts.std = 0; % Temperature - Standard deviation of distribution - replace with other parameters for method of moments
    PBM_val.Xs.mu = 0; % Reaction extent
    PBM_val.Xs.std = 0; % Reaction extent - 

    


%% Create a mesh


    %HX Section
    zs = linspace(0,reactor.H,Nz+1);
    rs = linspace(-reactor.R, reactor.R, Nr+1);

    [rr, zz] = meshgrid(rs, zs);
    zzz = ones(size(zz));
    meshopts.shape = 'Cylinder';
    meshopts.bottom_bound_row = true; %Create a row below mesh for initial value
    meshopts.boundary_names = {'Inlet', 'Outlet', 'Wall_left', 'Wall_right'}; %For plotting
    meshopts.boundary_curves = {[], [], [], []}; 

    mesh = RegularMesh(rr, zz, zzz, meshopts);
    

    %Transition




    %Reaction zone


    %Plot mesh
    
%% Setup plots

    %Plot inlet bubble size distribution
    figure();

    plot(100 .* mms, Fb_i, 'k-', 'LineWidth', lw); 
    xlabel('Bubble Radius (cm)'); ylabel('Rel. Number');
    grid on; grid minor; axis square;
    title('Initial Bubble Size Distribution (BSD)');

%% Calculate finite difference coefficients

    %Generate Mesh
    mesh = InitialValue(mesh, PBM_val);
    
    %Determine spatial finite difference stencil and coefficients
    fdopts.order = 2;
    fdopts.direct_bias = 0;
    mesh = FDCoeffs(mesh, 1, fdopts);

    %Determine bubble size finite different stencil and coefficients
    mfdopts.order = 2;
    mfdopts.direct_bias = 0;
    mmesh.type = '1D';
    mmesh.xsc = mms;
    mmesh.xsb = mbs;
    mmesh = FDCoeffs(mmesh, 1, mfdopts);

    %Model settings
    params.debug = true;
    params.coalescence = true;
    params.breakage = false;
    params.heat = true;
    params.reaction = false;

    %Define physical constants
    params.R = R;
    params.g = g;
    params.M_CH4 = 0.016014;
    params.M_H2 = 0.002016;
    params.M_Ar = 0.039948;

    %Define parameters
    params.scheme = 'MUSCL';
    params.mu_art = 1E-3;
    params.uz = 0.3;
    params.dz = reactor.H/Nz;
    params.h = reactor.H;
    params.mesh = mesh;
    params.rmesh = mmesh;
    params.N_cells = mesh.N_cells;
    params.Nz = Nz;
    params.bottom_inds = find(mesh.volcell_cents(:,2) < min(mesh.volcell_cents(:,2)) + 1E-6); %move to params
    params.top_inds = find(mesh.volcell_cents(:,2) > max(mesh.volcell_cents(:,2)) - 1E-6);
    params.N_rs = Nrs;
    params.rbs = rbs;
    params.rms = rms;
    params.Nms = Nms;
    params.mbs = mbs;
    params.mms = mms;
    params.alpha_g_max = 0.8;
    params.rho_L = rho_L;
    params.R = R;
    params.g = g;

    %Define heat transfer parameters
    params.T_Lz =  @(z) 1573.15 .* ones(size(z));
    

    %Calculate critical diameters for each bubble size for heat transfer
    Ts = linspace(reactor.T_min, reactor.T);
    ps = p_orifice .* ones(size(Ts));
    alphas = methane.therm_cond(Ts)./(methane.density(ps, Ts) .* methane.Cp(Ts)./methane.mol_mass);
    alpha_min = min(alphas);
    Fo_crit = 2;
    Lcs = params.rms;
    ubs = 
    t_ht = (Fo_crit .* Lcs.^2)./alpha_min;
    t_adv = params.h./
    
    %Initialize fields
    params.p_z =  InitializePressure(mesh, params);



    %Calculate volume of each spatial cell
    params.Vsz = pi .* reactor.R^2 .* (mesh.yy(2:end, :) - mesh.yy(1:end-1, :)); %m3 - volume of z cells

    %Add material properties
    params.tin = tin;
    params.methane = methane;
    params.hydrogen = hydrogen;


    %Determine which combinations of bubbles can coalesce
    params = IdentifyCoalescencePartners(params);
    
    %Generate equation matrix
    Mat_opts.z_dir = 2; % input number to meshgrid corresponding to z-direction (gravity) - 1 = first (columns), 2 = second (rows), 3 = third (depth)
    params.FD_mat = PBM_Mat(mesh, Mat_opts);

    %Create initial "volumes"
    N_volumes = mesh.N_cells .* Nms;
    params.cellinds = repmat([1:mesh.N_cells], Nms, 1); params.cellinds = params.cellinds(:); %Specifies which spatial cell the volume maps to
    params.xinds = repmat([1:Nms]', mesh.N_cells, 1);
    params.zinds = repmat([1:Nz], Nms, 1); params.zinds = params.zinds(:);

    %Define updated call spatial vectors
    cents_x = repmat(mesh.volcell_cents(:,1)', Nms, 1); cents_x = cents_x(:);
    cents_y = repmat(mesh.volcell_cents(:,2)', Nms, 1); cents_y = cents_y(:);
    params.cents_x = cents_x;
    params.cents_y = cents_y;

    %Calculate cell values
    params.p_orifice = p_surf + rho_L*g*reactor.H;
    params.nRTs = params.p_orifice .* (4/3) .* pi .* rms.^3; %Constant
    params.p_func = @(z) p_surf + rho_L.*g.*(H-z);
    params.ps = p_surf + rho_L*g*(reactor.H-params.cents_y);
    params.Ts = reactor.T .* ones(size(params.cents_y));
    params.Xs = (params.cents_y/reactor.H);
    params.ms = repmat(params.mms, 1, params.Nz); params.ms = params.ms';
    params.ns_i = params.ms./params.M_CH4;
    params.rhos_l = tin.density(params.Ts);
    params.sigmas = tin.surf_tension(params.Ts);

    %Calculate bubble velocity lookup table
    params.m_min = min(params.ms(:)); params.m_max = max(params.ms(:));
    params.p_min = min(params.ps(:)); params.p_max = max(params.ps(:));
    params.T_min = min(params.Ts(:)); params.T_max = max(params.Ts(:));
    params.V_max = (params.m_max .* params.R .* params.T_max)./(params.p_min .* params.M_H2);
    params.V_min = (params.m_min .* params.R .* params.T_min)./(params.p_max .* params.M_CH4);
    params.Ds_ubs_bounds = [(6 .* params.V_min/pi).^(1/3), (6 .* params.V_max/pi).^(1/3)];
    params.Ts_ubs = 1000:20:1400 + 273.15;
    params.Ds_ubs = linspace(params.Ds_ubs_bounds(1), params.Ds_ubs_bounds(2), 2.*params.Nms);
    params.ubs_lookup = CalcVelocities(params.Ds_ubs, params.Ts_ubs, params.tin, params.p_orifice);
    [params.Ts_ubs, params.Ds_ubs] = meshgrid(params.Ts_ubs, params.Ds_ubs);
    params.ub_func = @(Ts, Ds) interp2(params.Ts_ubs, params.Ds_ubs, params.ubs_lookup, Ts, Ds);

    %Calculate local gas superficial velocity
    
    params.uspfs = 1;
    

    %Create initial "volumes" and specify boundary conditions
    y0 = zeros(N_volumes, 1);
    bottom_inds = params.cents_y < min(mesh.volcell_cents(:,2)) + 1E-6;
    %y0(bottom_inds) = 1;
    Fi_bottom = repmat(Fb_i, 1, Nr);
    y0(bottom_inds) = Fi_bottom;


    %Solve ODE
    t_start = cputime;
    tspan = [0, t_end];
    tspan = 0:0.01:t_end;
    odeopts = odeset('InitialStep', 1e-6, 'MaxStep', 1e-1, 'RelTol', 1e-6, 'AbsTol', 1e-6);
    %[T,Y] = ode89(@(t,y) PBM_ode(t,y,params), tspan, y0);
    [T,Y] = ode45(@(t,y) PBM_ode(t,y,params), tspan, y0);
    t_end = cputime;
    t_run = abs(t_end - t_start);
    fprintf('ODE Solution Complete (t = %0.2f s)\n', t_run);

    %odeopts = odeset('Jacobian', @(t,y) PBM_jac(t,y, params));
    %[T,Y] = ode15s(@(t,y) PBM_ode(t,y,params), tspan, y0, odeopts);


    %Plot results
    

    


    

    %Interpolate result to vertices
    if Nr > 1

        %Show final value
        xind = 3;
        inds = xind:params.N_rs:((params.N_cells-1) * params.N_rs + xind);
        Yf = Y(end,inds);
        Yf = reshape(Yf, Nz, Nr);
    
        if Nr == 1 %1D
            rsc = 1;
        else %2D
            rsc = mesh.volcell_cents(:,1);
            rsc = reshape(rsc, Nz, Nr);
        end
        zsc = mesh.volcell_cents(:,2);
        zsc = reshape(zsc, Nz, Nr);
        rs = unique(round(rsc(:),4));
        zs = unique(round(zsc(:),4));
        [rsc,zsc] = meshgrid(rs, zs);
        Yfi = interp2(rsc, zsc, Yf, rr, zz, 'spline');

        figure()
        contourf(rr, zz, Yfi);
        colorbar(); clim([0,1])
        axis equal;
    else

        %Plot front
        figure();
        xind = 3;
        inds = xind:params.Nms:((params.N_cells-1) * params.Nms + xind);
        zsc = mesh.volcell_cents(:,2);
        zsc = reshape(zsc, Nz, Nr);
        subplot(1,2,1);
        Yf = Y(end,inds);
        plot(zsc, Yf)

        %Plot bubble size distribution
        subplot(1,2,2);
        
        Yf = Y(end,:);
        zsc = mesh.volcell_cents(:,2);
        zsc = repmat(zsc, 1, Nms);
        xsc = repmat(mmesh.xsc, mesh.N_cells,1);
        Yf = reshape(Yf,Nms, mesh.N_cells); Yf = Yf';

        surf(100.*zsc, 1000.*xsc, Yf);
        ylabel('Bubble Mass (\mug)')
        xlabel('Vertical Position (cm)')
        zlabel('Numerical Density');
        axis square; grid on; grid minor;
        

        %Animation plot
        if animate
            figure('Units', 'Normalized', 'OuterPosition', [0, 0, 1, 1]);
            for it = 1:length(T)
    
                t = T(it);
                Yit = Y(it,:);
                sc = mesh.volcell_cents(:,2);
                zsc = mesh.volcell_cents(:,2);
                zsc = reshape(zsc, Nz, Nr);
                zsc = repmat(zsc, 1, Nms);
                xsc = repmat(mmesh.xsc, mesh.N_cells,1);
                Yit = reshape(Yit,Nms, mesh.N_cells); Yit = Yit';
    
                surf(100.*zsc, 1000.*xsc, Yit);
                title(sprintf('t = %0.4f s', t));
                ylabel('Bubble Mass (\mug)')
                xlabel('Vertical Position (cm)')
                zlabel('Numerical Density');
                axis square; grid on; grid minor;
                view([45, 45])

                %Save result
                if saveanimation

                    % Capture the plot as an image 
                    frame = getframe(gcf); 
                    im = frame2im(frame); 
                    [imind,cm] = rgb2ind(im,256); 
                    
                    % Write to the GIF File 
                    if it == 1 
                      imwrite(imind,cm,fname,'gif', 'Loopcount',inf, 'DelayTime',gifdelay); 
                    elseif mod(it,10) == 0
                      imwrite(imind,cm,fname,'gif','WriteMode','append', 'DelayTime',gifdelay); 
                    end 
                    
                end
    
                pause(0.01);


            end
        end

        %

        x=1;
    end


    

    


end


function mesh =InitialValue(mesh, init)
    for ic = 1:mesh.N_cells
        mesh.volcells{ic}.val = init;
    end

end


function params = IdentifyCoalescencePartners(params)

    %Allocate output matrix
    params.coalesce_partners = cell(params.Nms, 1);
    params.coalesce_etas = cell(params.Nms, 1);

    %Iterate through representative sizes
    for xind = 1:params.Nms

        %Pull brackets
        mi = params.mms(xind);
        mb_low = params.mbs(xind); %Lower bound for this bracket
        mb_hi = params.mbs(xind+1); %Upper bound for this bracket

        %Determine adjecent cells
        if xind == 1
            mi_low = params.mms(xind);
        else
            mi_low = params.mms(xind-1);            
        end

        if xind == params.Nms
            mi_hi = params.mms(xind);
        else
            mi_hi = params.mms(xind+1);
        end

        %Pull cell value

    
        %Iterate through all possible size combinations
        coalesce_partners = [];
        eta_ijk = [];
        bias = [];
        for j = 1:params.Nms
            for k = j:params.Nms
               
                %Determine 
                mmj = params.mms(j);
                mmk = params.mms(k);
    
                mjk = mmj + mmk; %mass if they coalesce
    
                %Determine if the coalesced bubble would fall in this
                %   bracket
                if mjk >= mb_low && mjk < mb_hi

                    if isempty(coalesce_partners)
                        coalesce_partners = [j,k];
                    else
                        coalesce_partners(end+1,:) = [j,k];
                    end

                    %Calculate distribution coefficient for pair
                    if mjk < mi
                        bias(end+1) = -1;
                        eta_ijk(end+1) = (mjk - mi_low)./(mi - mi_low);
                    elseif mjk >= mi
                        bias(end+1) = 1; %1 indicates upper neighbor shares it
                        eta_ijk(end+1) = (mi_hi - mjk)./(mi_hi - mi);
                    end

                end

                

            end
        end


        


        %Log result
        params.coalesce_partners{xind} = coalesce_partners; %Lists bubble size pairs that coalesce into a specific pivot
        params.coalesce_etas{xind} = eta_ijk; %Indicates fraction assigned to main pivot
        params.coalesce_bias{xind} = bias; %Indicates which neighbor the value is shared with - 1 = hi, -1 = low
    end

end



function ubs = CalcVelocities(Ds, Ts, tin, p_orifice)

    %Physical constants
    g = 9.81;
    R = 8.3145;
    M_CH4 = 0.016;

    %Iterate through all combinations
    for it = 1:length(Ts)

        %Pull temperature
        T = Ts(it);

        %Calculate properties
        rho_G = (p_orifice .* M_CH4)./(R .* T);
        rho_L = tin.density(T);
        mu_L = tin.dyn_visc(T);
        sigma = tin.surf_tension(T);


        %Iterate through diameters
        for id = 1:length(Ds)

            %
            Db = Ds(id);

            %Calculate eotvos number
            Eob = ((rho_L - rho_G) .* g .* Db.^2)./sigma;
            
            %Calculate velocity
            ub = CalcVelocity(Db, rho_L, mu_L, sigma, Eob);

            %Log result
            ubs(id, it) = ub;

        end
    end


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



%Resolve pressure field in the domain
% - Account for non-uniform temperature and densities
function p_z = InitializePressure(mesh, params)

    
    


end




function temps = InitializeTemperature()


end



%A function to 
function outputs = TemperatureCharacteristics()

    %Iterate through each mass


        %Iterate through each temperature pivot



end


function dydt = TempOde(t,y,params)

    %Pull differential variables
    z = y(1);
    T = y(2);

    %Determine
    T_L = params.T_L(z);




    p = params.p_z(z);



    %Derivatives
    dTdt = 0;
    dzdt = 0;


    %Define output vector
    dydt = [dzdt; dTdt];


end