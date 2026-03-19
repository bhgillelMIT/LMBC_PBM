function results = PBM_v3(inputs)

%% Setup

    %Initialization Message
    fprintf('PBM - Initializing...\n\n')

    %Make parameters global
    global params

    %Clean up
    close all

    %Debug
    debug = false;



    %Plot settings
    lw = 2;
    fs = 18;

    %Setup folders
    addpath('src/');
    folders = SetupFolders();    

    %Material Properties
    addpath('Material Properties/');
    load tin.mat
    load hydrogen.mat
    load carbon.mat
    load methane.mat
    load argon.mat

    %Update values
    tin.dyn_visc = @(T) 0.00031 .* exp(6.171./(0.008314 .* T));

    %Add paths
    addpath('src/coalesce/');
    addpath('src/break/');
    addpath('Data/');
    addpath('Data/Breakage/');
    addpath('odesrc/');

    %Physical constants
    R = 8.3145;
    g = 9.81;

    



%% Inputs

    %Update user
    fprintf('-- Loading Inputs.\n')

    %Load inputs
    if nargin < 1
        inputs = PBM_inputs();
    end

    %Extract inputs
    animate = inputs.animate;
    fsolve_opts = inputs.fsolve_opts;
    mmesh = inputs.mmesh;
    reactor = inputs.reactor;
    disc = inputs.disc;
    sim = inputs.sim;
    sol = inputs.sol;
    inlet = inputs.inlet;
    initial = inputs.initial;
    liquid = reactor.liquid;

    %Calculate reactor flow rate
    reactor.V_dot = reactor.u_spf_orifice .* reactor.Ac; %m3/s
    reactor.V_dot_Lmin = 1000 .* reactor.V_dot .* 60;

    %Handle singler layer case
    [reactor, disc, sol] = PBM_singlelayer(reactor, disc, sol); %Adjust properties for case of a single cell

    %Parse inputs
    
    

%% Define mass mesh and boundary conditions

    %Update usr
    fprintf('-- Create mass mesh (Nm = %d). \n', disc.Nms)

    %Calculate gas density at the orifice
    reactor.T_bar = (reactor.T + reactor.T_min)/2;
    reactor.rho_L_bar = reactor.liquid.density(reactor.T_bar);
    p_orifice = reactor.p_surf + reactor.rho_L_bar * g * reactor.H;
    T_orifice = reactor.T_orifice;
    rho_gas_orifice = (p_orifice * reactor.M_gas)./(R * reactor.T_orifice);

    %Generate mass mesh (mmesh)
    [disc, mmesh] = PBM_mmesh(disc, mmesh, reactor, p_orifice, rho_gas_orifice, debug);

    %Define PBM data storage array
    PBM_val.eg = 0; % Gas fraction
    PBM_val.uz_L = 0; % Liquid vertical velocity
    PBM_val.Fbs = zeros(1,disc.Nrs); % Bubble volume - Discrete quantities
    PBM_val.Ts.mu = initial.T; % Temperature - Mean
    PBM_val.Ts.std = 0; % Temperature - Standard deviation of distribution - replace with other parameters for method of moments
    PBM_val.Xs.mu = 0; % Reaction extent - 
    PBM_val.Xs.std = 0; % Reaction extent - 

%% Create a mesh

    %Update user
    fprintf('-- Creating spatial mesh (Nz = %d).', disc.Nz)

    %HX Section
    zs = linspace(0,reactor.H,disc.Nz+1);
    rs = linspace(-reactor.R, reactor.R, disc.Nr+1);

    [rr, zz] = meshgrid(rs, zs);
    zzz = ones(size(zz));
    meshopts.shape = 'Cylinder';
    meshopts.bottom_bound_row = true; %Create a row below mesh for initial value
    meshopts.boundary_names = {'Inlet', 'Outlet', 'Wall_left', 'Wall_right'}; %For plotting
    meshopts.boundary_curves = {[], [], [], []}; 

    mesh = RegularMesh(rr, zz, zzz, meshopts);
    

%% Define Boundary Conditions

    %Update usr
    fprintf('-- Boundary Conditions (u_spf = %0.3f m/s). \n', reactor.u_spf_orifice)

    %Generate BCs
    [Nb_i, Fb_i] = PBM_BCs(reactor, mesh, inlet,  disc, zs, mmesh.rms, mmesh.Vms, T_orifice, reactor.liquid, p_orifice,  inputs, fsolve_opts);

    %Plot inlet bubble size distribution
    if debug
        figure();
    
        plot(100 .* mmesh.rms, Fb_i, 'k-', 'LineWidth', lw); 
        xlabel('Bubble Radius (cm)'); ylabel('Num. density (1/(kg*m^3))');
        grid on; grid minor; axis square;
        title('Initial Bubble Size Distribution (BSD)');
    end

    

%% Generate Mesh

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
    mmesh.xsc = mmesh.mms;
    mmesh.xsb = mmesh.mbs;
    mmesh = FDCoeffs(mmesh, 1, mfdopts);

%% Model Inputs

    

    %Model settings
    params.debug = true;
    params.logresults = true;
    params.logperiod = 0.01; %s - space between recording 
    params.sol.orifice_BC_type = sol.orifice_BC_type; %1 = Dirichlet; 2 = Neumann
    params.sol.rel_tol = 1E-4;
    params.sol.abs_tol = 1E-4;
    params.dt = disc.dt;
    params.dt_max = 0.01;
    params.dt_min = 0.000001; 
    params.t_f = sim.t_end;
    params.fsolve_opts = fsolve_opts; 
    params.output.animate = animate.active;
    params.output.saveanimation = animate.save;
    params.output.folder = folders.sol;
    params.output.gifdelay = inputs.animate.gifdelay;
    params.folders = folders;
    params.disc = disc;

    %Solution parameters
    params.sol.debug = true;  
    params.sol.type = sol.type; %Type of solver to use: 'direct' solves all equations simultaneously; 'segregated' use strang operator splitting to separate advection from source terms
    params.sol.scheme = sol.scheme ; %'Upwind'; %Advection scheme - options: 'MUSCL'; 'Upwind'
    params.sol.mu_art = 1E-6;
    params.sol.src_delay = sol.src_delay; %second - time delay for the start of source terms
    params.sol.sep_layer = true; %boolean - determines if source terms are solved all at once (false), or separately (true).
    params.sol.solve_ub = sol.solve_ub;
    params.sol.ub_manual = sol.ub_manual;
    params.sol.solve_details = true;
    params.sol.solve_detail_its = sol.solve_detail_its;
    params.sol.single_layer = sol.single_layer;

    %Characteristic parameters
    params.chars.debug = false;
    params.chars.fig_dir = 'Figures/Characteristics/';
    params.chars.data_dir = folders.char;
    params.chars.print = true;
    params.chars.type = 'LC'; % 'LC' = lumped capacitance (simple); 'SBM' = Single Bubble Model'; '1D' = 1D spherical shell resistance model using LC = R/3
    params.chars.LC_mode = 'Cond'; % 'Conv' = low biot number; 'Cond' = high biot
    params.chars.load = true;
    params.chars.interp_method = 'linear';
    params.chars.N_Xis = 4;
    params.chars.dT_cutoff = 10; %K - threshold for merging temperature bins
    params.chars.dX_cutoff = 0.02; % percent - threshold for merging conversion bins
    params.chars.unique_thresh = 0.0001;

    %Source term parameters
    params.src.debug = true;
    params.src.solve_eps = inputs.src.solve_eps; %Boolean - toggle between solving the turbulent kinetic energy dissipation rate (true) or not (false)
    params.src.eps_manual = inputs.src.eps_manual; %
    params.src.solve_alphag = inputs.src.solve_alphag;
    params.src.alphag_manual = inputs.src.alphag_manual;
    params.src.its = 0; %Tracker for number of calls to the source functions (coalescence and breakage)
    

    %Coalescence parameters
    params.coalesce.active = inputs.src.coalesce_active;
    %params.coalesce.constant = true; %
    params.coalesce.eddy = true; %Eddy coalescence - true = enabled; false = disabled
    params.coalesce.wake = false; %Wake entrainment coalescence - true = enabled; false = disabled
    params.coalesce.rise = true; %Coalescence from different rise velocities - true = enabled; false = disabled;
    params.coalesce.damper = 1; %Coalescence rate damper (to help with visualization) - SHOULD BE 1 OTHERWISE
    params.coalesce.We_crit = 100; %Critical weber number for bubble breakup - from Wang et al. 2005
    params.coalesce.constant_freq = 1; %Collisions per second
    params.coalesce.m_src = zeros(1, 10000);
    params.coalesce.m_snk = zeros(1, 10000);
    params.coalesce.model = inputs.src.coalesce_model;
    params.coalesce.constant_rate = inputs.src.coalesce_constant_rate;

    %Breakage parameters
    params.break.active = inputs.src.breakage_active;
    params.break.debug = false;
    params.break.eddy = true;
    params.break.surf = true;
    params.break.damper = 1; %Breakage rate damper to help with visualization - SHOULD BE 1 OTHERWISE
    params.break.b_star = 100; %Model parameter - Wang et al 2005
    params.break.m_star = 6.0; %Model parameter - 
    params.break.folder = folders.break;
    params.break.loadfile = true; %Determines if a new interpolation file is first generated - set to true to use the existing file referenced 
    params.break.model = inputs.src.breakage_model;
    params.break.type = 'non-uniform'; %Type of breakage - 'uniform' = bubbles are uniformly distributed among smaller sizes; 'non-uniform' = the an exact bubble size distribution, this is slower.
    params.break.int_method = 'trapz'; %Type of integration to use for evaluating integrals in breakage kernels - Options: 'gauss' for gaussian quadrature, 'trapz' for trapezoidal integration
    params.break.badd = NaN;
    params.break.bsub = NaN;
    params.break.m_src = zeros(1, 10000);
    params.break.m_snk = zeros(1, 10000);
    params.break.model = inputs.src.breakage_model;
    params.break.constant_rate = inputs.src.breakage_constant_rate;


    %Heat transfer parameters
    params.heat.active =  sol.heat.active;
    params.heat.source_mode = 'simple'; % Options: 'simple' or 'adaptive'; Details: 'simple' allocation proportional to original allocation for the sink representative mass; 



    %Reaction parameters
    params.react.active = sol.react.active;

    

    %Define physical constants and materials
    params.T0 = 298.15;
    params.p0 = 101325;
    params.R = R;
    params.g = g;
    params.liquid = liquid;
    params.gas = 1; %IMPLEMENT
    params.methane = methane;
    params.hydrogen = hydrogen;
    params.carbon = carbon;
    %params.tin = tin;
    params.argon = argon;
    params.M.CH4 = 0.016014;
    params.M.H2 = 0.002016;
    params.M.Ar = 0.039948;   

    
    



    %Define parameters
    params.zms = (mesh.yy(1:end-1,1) + mesh.yy(2:end,1))./2;
    params.dz = reactor.H/disc.Nz;
    params.reactor = reactor;
    params.mesh = mesh;
    params.rmesh = mmesh;
    params.mmesh = mmesh;
    params.disc = disc;
    params.N_cells = mesh.N_cells;
    params.Nr = disc.Nr;
    params.Nz = disc.Nz;
    params.bottom_inds = find(mesh.volcell_cents(:,2) < min(mesh.volcell_cents(:,2)) + 1E-6); %move to params
    params.top_inds = find(mesh.volcell_cents(:,2) > max(mesh.volcell_cents(:,2)) - 1E-6);
    %params.N_rs = Nrs;
    params.rbs = mmesh.rbs;
    params.rms = mmesh.rms;
    params.dbs = 2 .* params.rbs;
    params.dms = 2 .*  params.rms;
    params.dbds = diff(params.dbs);
    params.Vms = mmesh.Vms; %(4/3) .* pi .* params.rms.^3;
    params.Nms = disc.Nms;
    params.mbs = mmesh.mbs;
    params.mds = diff(mmesh.mbs); %spacing of mass brackets
    params.mms = mmesh.mms;
    params.nms = mmesh.nms; %mols - number of moles originally
    params.nbs = mmesh.nbs;
    params.alpha_g = zeros(1, params.Nz);
    params.Vcells = reactor.Ac .* diff(params.mesh.yy(:,1));
    params.p_surf = reactor.p_surf;
    params.alpha_g_max = 0.8;
    params.X_Ar = 0;
    params.M_gas_i = params.argon.mol_mass * params.X_Ar + params.methane.mol_mass * (1-params.X_Ar);
    %params.N_volumes = params.N_cells * params.Nms;
    
    %Define heat transfer parameters
    params.dT_end_max = 20; %K - maximum spacing of points at the 
    params.h_conv_const = true;
    params.T_Lz =  @(z) reactor.T .* ones(size(z)); %K - Liquid temperature profile function
    params.rho_L_T = @(T) reactor.liquid.density(T); %@(T) 6979 - 0.652 * (T - 505.08); %kg/m3 - Liquid density as a function of temperature
    params.T_min = reactor.T_min; %(T_mu - 4*T_std);
    params.T_mu_i = reactor.T_gas_i_mu;
    params.T_std_i = reactor.T_gas_i_std;
    params.T_liq = reactor.T;
    params.M_gas = params.M_gas_i ; %methane.mol_mass;
    params.h_conv = @(T) 5 .* ones(size(T));
    params.Cp = @(T) 1000 .* ones(size(T));

    %Determine if system is isothermal
    T_diffs = params.T_Lz(params.zms) - params.T_Lz(params.zms(1));
    params.isothermal = all(T_diffs < 0.1);

    %Define and initialize fluid parameter
    params.turb.eps = zeros(size(params.zms));
    params.turb.k = zeros(size(params.zms));
    params.mus = liquid.dyn_visc(params.T_Lz(params.zms));
    params.rhos = liquid.density(params.T_Lz(params.zms));
    params.sigmas = liquid.surf_tension(params.T_Lz(params.zms));
    params.nus = params.mus(:)./params.rhos(:);

    %Calculate critical diameters for each bubble size for heat transfer
    Ts = linspace(reactor.T_min, reactor.T);
    ps = p_orifice .* ones(size(Ts));
    alphas = methane.therm_cond(Ts)./(methane.density(ps, Ts) .* methane.Cp(Ts)./methane.mol_mass);
    alpha_min = min(alphas);
    Fo_crit = 1;
    Lcs = params.rms;
    T_bar = params.T_Lz(params.reactor.H/2);
    ubs = CalcVelocities(2.*Lcs, T_bar, liquid, p_orifice, fsolve_opts);
    t_ht = (Fo_crit .* Lcs.^2)./alpha_min;
    t_adv = params.reactor.H./ubs;

    
    %Initialize fields
    params.p_z =  InitializePressure(mesh, params);
    params.p_orifice = params.p_z(0);

    %Calculate volume of each spatial cell
    params.Vsz = pi .* reactor.R^2 .* (mesh.yy(2:end, :) - mesh.yy(1:end-1, :)); %m3 - volume of z cells


    %Determine which combinations of bubbles can coalesce
    params = IdentifyCoalescencePartners(params);
    
    %Generate equation matrix
    Mat_opts.z_dir = 2; % input number to meshgrid corresponding to z-direction (gravity) - 1 = first (columns), 2 = second (rows), 3 = third (depth)
    params.FD_mat = PBM_Mat(mesh, Mat_opts);

    %Define updated call spatial vectors
    cents_x = repmat(mesh.volcell_cents(:,1)', disc.Nms, 1); cents_x = cents_x(:);
    cents_y = repmat(mesh.volcell_cents(:,2)', disc.Nms, 1); cents_y = cents_y(:);
    params.cents_x = cents_x;
    params.cents_y = cents_y;

    
    
    %Load or resolve temperature characteristics
    params.LC_trail_type = 'Initial'; %How trailing points are designated: 'Final' = the final temperature is defined to maintain a known resolution at the end state; 'Initial' = initial temperature defined as somewhere on the slowest line, then final temperature resolved
    params.LC_trail_spacing = 'geo'; %How the trailing points are spaced: 'lin' = equal spacing; 'geo' = geometric spacing; 'log' = logarithmic spacing
    params.LC_trail_geo_r = 0.7;
    params.SBM_folder = 'SBM Characteristics/Demo6/';
    params = LoadTempChars(params); %Load or resolve T/X characteristics
    

    %Create initial "volumes"
    [y0, params, params.N_volumes, Fb_i] = InitializeVolumes(params, mesh, disc, Fb_i);
    %N_volumes = mesh.N_cells .* disc.Nms;
    % params.cellinds = repmat([1:mesh.N_cells], disc.Nms, 1); params.cellinds = params.cellinds(:); %Specifies which spatial cell the volume maps to
    % params.xinds = repmat([1:disc.Nms]', mesh.N_cells, 1);
    % params.zinds = repmat([1:disc.Nz], disc.Nms, 1); params.zinds = params.zinds(:);

    %Calculate cell values
    params.p_orifice = reactor.p_surf + reactor.rho_L_bar*g*reactor.H;
    params.nRTs = params.p_orifice .* (4/3) .* pi .* mmesh.rms.^3; %Constant
    params.p_func = @(z) reactor.p_surf + reactor.rho_L_bar.*g.*(reactor.H-z);
    params.ps = reactor.p_surf + reactor.rho_L_bar*g*(reactor.H-params.zcs_rep);
    params.Ts = params.T_Lz(params.zcs_rep);%   reactor.T .* ones(size(params.cents_y));
    params.Tsz = params.T_Lz(unique(params.cents_y));
    params.Xs = (params.zcs_rep/reactor.H);
    params.ms = params.mms_rep; %repmat(params.mms, 1, params.Nz); params.ms = params.ms';
    params.ns_i = params.ms./params.M.CH4;
    params.rhos_l = liquid.density(params.Ts);



    %Calculate bubble velocity lookup table - only for non-reacting,
    %   isothermal case 
    params.ubs.m_min = min(params.ms(:)); params.ubs.m_max = max(params.ms(:));
    params.ubs.p_min = min(params.ps(:)); params.ubs.p_max = max(params.ps(:));
    params.ubs.T_min = min(params.Ts(:)); params.ubs.T_max = max(params.Ts(:));
    params.ubs.V_max = (params.ubs.m_max .* params.R .* params.ubs.T_max)./(params.ubs.p_min .* params.M.H2);
    params.ubs.V_min = (params.ubs.m_min .* params.R .* params.ubs.T_min)./(params.ubs.p_max .* params.M.CH4);
    params.ubs.Ds_ubs_bounds = [(6 .* params.ubs.V_min/pi).^(1/3), (6 .* params.ubs.V_max/pi).^(1/3)];
    params.ubs.Ts_ubs = 1000:20:1400 + 273.15;
    params.ubs.Ds_ubs = linspace(params.ubs.Ds_ubs_bounds(1), params.ubs.Ds_ubs_bounds(2), 2.*params.Nms);
    params.ubs.funcs = BubbleVelocities(params);
    

    %Calculate range of superficial velocities
    params.n_dot_orifice = (reactor.V_dot .* params.p_orifice)./(R .* reactor.T_gas_i_mu);
    params.u_spf_orifice = reactor.u_spf_orifice;
    params.V_dot_orifice = reactor.V_dot;
    params.V_dot_surf = (2 .* params.n_dot_orifice .* R .* reactor.T)./(params.p_surf);
    params.u_spf_surface = params.V_dot_surf/reactor.Ac;
    params.u_spfs = [reactor.u_spf_orifice, params.u_spf_surface];


    %Precalculate mean conversions, temperatures, and diameters 
    params.X_mu = 0.1 .* ones(params.Nz, params.Nms);
    %params.T_mus = ones(params.Nz, params.Nms);
    params.T_mu = params.T_mu_i .* ones(params.Nz, params.Nms); %reshape(params.Ts, params.Nz, params.Nms);
    params.d_mu = zeros(size(params.T_mu));
    params.uz_mu = zeros(size(params.T_mu));
    params.V_mu = zeros(size(params.T_mu));
    params.Vb_mu = zeros(size(params.rbs));
    params.n_mu = params.V_mu;
    params.nb_mu = params.Vb_mu;
    params.db_mu = params.Vb_mu;
    params.ug = zeros(params.Nz, 1);
    params.Ns_z_zero = zeros(1, params.Nz); %1/m3 - Total Numeric density for each spatial cell
    params.Ns_m_zero = zeros(1, params.Nz * params.Nms); %1/m3 - Total Numeric density for each mass cell
    params.Ns_T_zero = zeros(1, length(find(params.Xinds == 1))); %1/m3 - Total Numeric density for each temperature cell
    params.h = zeros(params.N_volumes, 1); %m3/s - source term storage vector - preallocating 
    params.zinds_m = repmat([1:params.Nz], params.Nms, 1); params.zinds_m = params.zinds_m(:);
    params.xinds_m = repmat([1:1:params.Nms], 1, params.Nz); params.xinds_m = params.xinds_m(:);

    %Breakage parameters
    params.break.delta = 0.01;
    params.break.N_lambdas = 300; %Number of eddy diameters to use for integration in eddy breakage calculations                                                                         
    params.break.dfv_surf = 0.025; %Resolution of breakup fraction (fv) at larger sizes
    params.break.dfv_mid = 0.01; %Resolution of breakup fration for intermediate sizes
    params.break.dfv_cap = 0.0025; %Resolution of breakup fraction for small sizes 
    params.break.fv_thresh_low = 0.05; %Threshold between small sizes and intermediate sizes
    params.break.fv_thresh_hi = 0.15; %Threshold between intermediate and large sizes
    params.break.N_fv_min = 10; %Minimum number of breakup fractions to use for the capillary pressure limited region
    params.fvs_norm = [0:params.break.dfv_cap:params.break.fv_thresh_low,...
                       (params.break.fv_thresh_low+params.break.dfv_mid):params.break.dfv_mid:params.break.fv_thresh_hi,...
                       (params.break.fv_thresh_hi+params.break.dfv_surf):params.break.dfv_surf:0.5]; %Normal breakup fraction values to interpolate results onto and to do subsequent integration with
    params.fvs_norm_all = [params.fvs_norm, 1 - fliplr(params.fvs_norm(1:end-1))]; %Full spectrum
    params.Pbs_norm = zeros(params.break.N_lambdas, length(params.fvs_norm)); %
    params.break.lambda_rats = linspace(1E-6,1, params.break.N_lambdas);
    params.fvc_func = Calculatefvcfunc(params.break.lambda_rats);
    params.break.bfd_zero = zeros(size(params.fvs_norm));
    params.break.fvih_saveoutput = true;
    params.break.fvih_loadoutput = true;
    params.fvih_max_funcs = Loadfvihfuncs(params, disc);

    %Create Bubble Size Distribution (BSD) interpolation function   
    params.break.interp = true;
    params.break.eps_range = [0, 5];
    params.break.d_range = [0.001, 0.05];
    params.break.N_u_spfs = 25;
    params.break.N_d_interp = 2;
    params.break.N_eps_interp = 2;
   % params.break.N_lambdas = 50;
    params.break.N_eps = params.break.N_u_spfs;
    params.break.N_ds = 25; %1.*params.Nms;
    params.break.filename = sol.break_file;
    params = LoadBreakFile(params);

    %Debug analyis
    if params.debug
        %PBM_Analysis(params); %Generate plots for evaluating the model inputs
    end


    %Solve bin temperatures
    params = CalcCharTemps(y0, params); %Resolve the middle bin temperature based on those characteristics
 


    %Allocate storage for detailed calculations to enable steps  without
    %solutions

    %Calculate initial mass in the system
    %params.Vcells_rep = repmat(params.Vcells, 1, params.Nms); params.Vcells_rep = params.Vcells_rep';
    %params.Vcells_rep = params.Vcells_rep(:);
    %params.mms_rep = repmat(params.mms', params.Nz, 1);
    params.m_total = zeros(1, 10000);
    params.m_total(1) = sum(params.mms_rep .* y0 .* params.Vcells_rep);

    
    

    %Define boundary condition
    params.N_dot_o = Nb_i/reactor.Ac;
    params.N_dot_os = CalcInletFluxes(params); %ISSUE TO RESOLVE

    %Debug tests
    %BreakageTest(params);
    %SourceTest(params);



%% Run Model

    %Debug output
    if params.debug
        fprintf('PBM Start (t = 0.0000000 s)')
    end    

    %Solve ODE
    params.its = 1;
    t_start = cputime;
    tspan = 0:0.01:sim.t_end; %tspan = [0, t_end];
    odeopts = odeset('InitialStep', 1e-3, 'MaxStep', params.dt_max, 'RelTol', 1e-5, 'AbsTol', 1e-5, 'Stats','on', 'OutputFcn', @PBM_output); %, 'NonNegative', 1:length(y0));    
    [T,Y] = PBM_solver(tspan, y0, params); % [T,Y] = ode45(@(t,y) PBM_ode(t,y,params), tspan, y0, odeopts);
    t_end = cputime; 
    t_run = abs(t_end - t_start); 
    fprintf('ODE Solution Complete (t = %0.2f s)\n', t_run);


    %Plot results
    %PBM_plot(Y, T, params)

    %Save results
    PBM_saveresults(T, Y, params)
    

    %Post process results
    results = PBM_postprocess(T, Y, params);


    

    % %Interpolate result to vertices
    % if disc.Nr > 1
    % 
    %     %Show final value
    %     xind = 3;
    %     inds = xind:params.N_rs:((params.N_cells-1) * params.N_rs + xind);
    %     Yf = Y(end,inds);
    %     Yf = reshape(Yf, Nz, Nr);
    % 
    %     if Nr == 1 %1D
    %         rsc = 1;
    %     else %2D
    %         rsc = mesh.volcell_cents(:,1);
    %         rsc = reshape(rsc, Nz, Nr);
    %     end
    %     zsc = mesh.volcell_cents(:,2);
    %     zsc = reshape(zsc, Nz, Nr);
    %     rs = unique(round(rsc(:),4));
    %     zs = unique(round(zsc(:),4));
    %     [rsc,zsc] = meshgrid(rs, zs);
    %     Yfi = interp2(rsc, zsc, Yf, rr, zz, 'spline');
    % 
    %     figure()
    %     contourf(rr, zz, Yfi);
    %     colorbar(); clim([0,1])
    %     axis equal;
    % else
    % 
    %     %Plot front
    %     figure();
    %     xind = 3;
    %     inds = xind:params.Nms:((params.N_cells-1) * params.Nms + xind);
    %     zsc = mesh.volcell_cents(:,2);
    %     zsc = reshape(zsc, disc.Nz, disc.Nr);
    %     subplot(1,2,1);
    %     Yf = Y(end,inds);
    %     plot(zsc, Yf)
    % 
    %     %Plot bubble size distribution
    %     subplot(1,2,2);
    % 
    %     Yf = Y(end,:);
    %     zsc = mesh.volcell_cents(:,2);
    %     zsc = repmat(zsc, 1, disc.Nms);
    %     xsc = repmat(mmesh.xsc, mesh.N_cells,1);
    %     Yf = reshape(Yf,disc.Nms, mesh.N_cells); Yf = Yf';
    % 
    %     surf(100.*zsc, 1000.*xsc, Yf);
    %     ylabel('Bubble Mass (\mug)')
    %     xlabel('Vertical Position (cm)')
    %     zlabel('Numerical Density');
    %     axis square; grid on; grid minor;
    % 
    % 
    %     %Animation plot
    %     if animate.active
    %         figure('Units', 'Normalized', 'OuterPosition', [0, 0, 1, 1]);
    %         for it = 1:length(T)
    % 
    %             t = T(it);
    %             Yit = Y(it,:);
    %             sc = mesh.volcell_cents(:,2);
    %             zsc = mesh.volcell_cents(:,2);
    %             zsc = reshape(zsc, disc.Nz, disc.Nr);
    %             zsc = repmat(zsc, 1, disc.Nms);
    %             xsc = repmat(mmesh.xsc, mesh.N_cells,1);
    %             Yit = reshape(Yit,disc.Nms, mesh.N_cells); Yit = Yit';
    % 
    %             surf(100.*zsc, 1000.*xsc, Yit);
    %             title(sprintf('t = %0.4f s', t));
    %             ylabel('Bubble Mass (\mug)')
    %             xlabel('Vertical Position (cm)')
    %             zlabel('Numerical Density');
    %             axis square; grid on; grid minor;
    %             view([45, 45])
    % 
    %             %Save result
    %             if params.output.saveanimation
    % 
    %                 % Capture the plot as an image 
    %                 frame = getframe(gcf); 
    %                 im = frame2im(frame); 
    %                 [imind,cm] = rgb2ind(im,256); 
    % 
    %                 % Write to the GIF File 
    %                 if it == 1 
    %                   imwrite(imind,cm,animate.fname,'gif', 'Loopcount',inf, 'DelayTime',animate.gifdelay); 
    %                 elseif mod(it,10) == 0
    %                   imwrite(imind,cm,animate.fname ,'gif','WriteMode','append', 'DelayTime',animate.gifdelay); 
    %                 end 
    % 
    %             end
    % 
    %             pause(0.01);
    % 
    % 
    %         end
    %     end
    % 
    %     %
    % 
    %     x=1;
    % end
    % 
    % 
    % 
    % 
    % 
    % 

    %Nested output function 

    % function status = PBM_output_nested(t, y, flag)
    % 
    %     %Allocate storage vectors
    %     persistent its 
    % 
    % 
    %     if isempty(its)
    %         its = 1;
    %     else
    %         its = its + 1;
    %     end
    % 
    %     if isempty(flag)
    %         flag = 'None';
    %     elseif strcmp(flag, 'done')
    %         its = 0; %Reset storage vector
    %     end
    % 
    %     %Print update
    %     fprintf('\nPBM (Iteration %d; flag = %s;) \n\n', its, flag)
    % 
    % 
    %     %Define status
    %     status = 0;
    % end
end


function mesh =InitialValue(mesh, init)
    for ic = 1:mesh.N_cells
        mesh.volcells{ic}.val = init;
    end

end


function params = IdentifyCoalescencePartners(params)

    fprintf('-- Identifying Coalescence Partners.\n');


    %Allocate output matrix
    params.coalesce_partners = cell(params.Nms, 1);
    params.coalesce_etas = cell(params.Nms, 1);
    params.coalesce_rats = cell(params.Nms, 1);

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

        %Identify max index to test 
        max_ind = xind; 

    
        %Iterate through all possible size combinations
        coalesce_partners = [];
        eta_ijk = [];
        bias = [];
        overrat = [];
        for j = 1:max_ind
            for k = j:max_ind
               
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
                    %Calculate distribution coefficient for pair
                    if mjk < mi
                        bias(end+1) = -1;
                        eta_ijk(end+1) = (mjk - mi_low)./(mi - mi_low);
                    elseif mjk >= mi
                        bias(end+1) = 1; %1 indicates upper neighbor shares it
                        if xind == params.Nms
                            eta_ijk(end+1) = 1;
                        else
                            eta_ijk(end+1) = (mi_hi - mjk)./(mi_hi - mi);
                        end
                    end




                    % if mjk < mi
                    %     bias(end+1) = -1;
                    %     eta_ijk(1, end+1) =  (mi - mjk)/(mi - mi_low); %(mjk - mi_low)./(mi - mi_low); %Smaller row is smaller cell (either mi_low, or mi depending on the location of mjk)
                    %     %eta_ijk(2, end) = (mjk - mi_low)./(mi - mi_low); 
                    % elseif mjk >= mi
                    %     bias(end+1) = 1; %1 indicates upper neighbor shares it
                    %     if xind == params.Nms
                    %         eta_ijk(1, end+1) = mjk/mi; %This is to preserve mass, not number. Older value of 1 
                    %         %eta_ijk(2, end) = 0; % Cannot allocate to bins that don't exist
                    %     else
                    %         eta_ijk(1, end+1) = (mi_hi - mjk)./(mi_hi - mi);
                    %         eta_ijk(2, end) = (mjk - mi)./(mi_hi - mi);
                    %     end
                    % end

                    %Sepcify ratios
                    overrat(end+1) = 1;

                end

                %Determine if they add to something bigger than the largest
                %bin
                if xind == params.Nms && mjk > mb_hi
                    coalesce_partners(end+1,:) = [j,k];
                    overrat(end+1) = mjk/params.mms(end);% Ratio over the mass of 
                    eta_ijk(1, end+1) = mjk/params.mms(end);
                    %eta_ijk(2, end) = 0;
                    bias(end+1) = 0;
                    x =1;
                % elseif xind == params.Nms && (mjk >= mb_low && mjk < mb_hi)
                %     overrat(end+1) = 1;
                end


                

            end
        end


        


        %Log result
        params.coalesce_partners{xind} = coalesce_partners; %Lists bubble size pairs that coalesce into a specific pivot
        params.coalesce_etas{xind} = eta_ijk; %Indicates fraction assigned to main pivot
        params.coalesce_bias{xind} = bias; %Indicates which neighbor the value is shared with - 1 = hi, -1 = low
        params.coalesce_rats{xind} = overrat;
        
    end

end






% function ub = CalcVelocity(Db, rho_L, mu_L, sigma, Eob)
% 
%     %Physical constants
%     g = 9.81;
% 
%     %Initial velocity guess
%     ub_guess = 0.3;
% 
%     ub_func = @(ub) sqrt((4*Db*g)./(3 * CalcCd(ub, Db, rho_L, mu_L, sigma, Eob))) - ub;
%     ub = fsolve(ub_func, ub_guess);
%     
%     
% 
% 
% end

% function Cd = CalcCd(ub, Db, rho_L, mu_L, sigma, Eob)
% 
%     %Calculate reynolds number
%     if ub > 0
%         Reb = (rho_L * ub * Db)/mu_L;
%     else
%         Reb = 1E-6;
%     end
% 
%     %Calculate drag coefficient
%     Cd_sph = (24/Reb)*(1 + 0.15*Reb^(0.687));
%     Cd_ell_cap = 8/3 * Eob/(Eob + 4);
%     Cd = max([Cd_sph, Cd_ell_cap]);
% 
% end



%Resolve pressure field in the domain
% - Account for non-uniform temperature and densities
function p_z = InitializePressure(mesh, params)

    %Calculate density at each height
    zs = linspace(0,params.reactor.H, 10 * params.Nz); zs = fliplr(zs);
    % zs_m = (zs(1:end-1) + zs(2:end))./2; zs_m = fliplr(zs_m);
    % dz = mode(diff(zs));
    % Ts_m = params.T_Lz(zs_m);
    % rhos_m = params.rho_L_T(Ts_m);
    
    %Integrate
    int_func = @(z) params.g .* params.rho_L_T(params.T_Lz(z));
    ps = zeros(1, length(zs));
    for iz = 1:length(zs)
        
        z = params.reactor.H - zs(iz);
        int_act = integral(int_func, 0, z);
        ps(iz) = params.p_surf + int_act;
 
    end

    %Define lookup table
    %ps = fliplr(ps);
    p_z = @(z) interp1(zs, ps, z);
    


    

end




function Tparams_output = InitializeTemperature(params)

    %Debug settings
    debug = true;

    %Define output vector 
    Tparams_output = cell(1, length(params.mms));

    %Pull relevant values
    Tparams.m_max = max(params.mms);

    %Iterate through size brackets and create initial points
    for im = 1:length(params.mms)

        %Pull values
        m = params.mms(im);
        r = params.rms(im);
        N_bs = 1000;

        %Initialize temp
        Tparams.h_conv_const = params.h_conv_const;
        Tparams.m = m;
        Tparams.im = im;
        Tparams.p = params.p_orifice;
        Tparams.p_surf = params.p_surf;
        Tparams.T_mu = params.T_mu_i;
        Tparams.T_std = params.T_std_i;
        Tparams.N_bs = N_bs;
        Tparams.T_liq = params.T_liq;
        Tparams.T_i = params.T_mu_i;
        Tparams.M_gas = params.M_gas;
        Tparams.h = params.h_conv;
        Tparams.Cp = params.Cp;
        Tparams.dt = params.dt;
        Tparams.t_f = params.t_f;
        Tparams.T_min = params.T_min;
        Tparams.core.N_pivots = 3; %Number of pivots per side
        Tparams.core.N_pivots_trail = 3;
        Tparams.core.N_stds = 3;
        Tparams.trail.dT = 25; %Minimum resolution to maintain
        Tparams.trail.dT_init = 5; %Initial spacing of trail points 
        Tparams.lead.dT = 75;
        Tparams = GenerateMovingPivots(Tparams, params);

        %Log distribution
        Tparams_output{im} = Tparams;

    end

    %Debug plot
    if debug

        pivots.Y = Tparams.Tps;
        pivots.Tms = (pivots.Y(1:Tparams.N_Ts-1) + pivots.Y(2:Tparams.N_Ts))./2;
        pivots.Fms = Tparams.Fps;
        plotconfig.lw = 2; plotconfig.fs = 18; plotconfig.ms = 18;

        PlotMovingPivots(Tparams, pivots, plotconfig); %Need to improve plot
    end

end

function fvc_func_poutput = Calculatefvcfunc(lambda_rats)
    
    fvcs = zeros(size(lambda_rats));
    for ir = 1:length(lambda_rats)
        lambda_rat = lambda_rats(ir);
        fvc_func = @(fvc) (1./216).*(lambda_rat).^9 - fvc .* (fvc.^(2/3) + (1-fvc).^(2/3) - 1).^3; %Use with fzero to find critical breakup fraction
        fvcs(ir) = fzero(fvc_func, [0,0.5]);
    end

    fvc_func_poutput = @(r) interp1(lambda_rats, fvcs, r);
    

end




% 
% function params = CalculateBSDs(params)
% 
%     %Calculate range of gas velocities
%     u_spf_min = params.u_spf_orifice;
%     u_spf_max = params.u_spf_surface;
%     u_spfs = linspace(0.9*u_spf_min, 1.1*u_spf_max, params.break.N_u_spfs);
% 
%     %Calculate range of bubble diameters
%     V_min = (params.nms(1) .* params.R .* params.T_mu_i)./params.p_orifice;
%     V_max = (2.*params.nms(end) .* params.R .* params.T_liq)./params.p_surf;
%     d_min = (6.*V_min/pi).^(1/3);
%     d_max = (6.*V_max/pi).^(1/3);
%     ds = logspace(log10(params.break.d_range(1)), log10(params.break.d_range(2)), params.break.N_ds);
% 
%     %Calculate range of turbulent kinetic energy dissipation (epsilon)
%     turb_epss = [0, logspace(log10(0.1), log10(params.break.eps_range(2)), params.break.N_u_spfs-1)]; %params.g * u_spfs';
% 
% 
% 
%     %Determine number of spatial cells to iterate through
%     if params.isothermal
%         Nz = 1;
%     else
%         Nz = params.Nz;
%     end
% 
%     %Iterate through spatial cells
%     %params.bs = cell(1, params.Nz);
%     params.beta_ratio = cell(1, params.Nz);
%     params.betas = cell(1, params.Nz);
%     for iz = 1:Nz
% 
%         %Create mesh grids
%         P = [2, 1, 3];
%         [ds_mg, turb_epss_mg, fvs_mg] = meshgrid(ds, turb_epss, params.fvs_norm_all);
%         ds_mg = permute(ds_mg, P); turb_epss_mg = permute(turb_epss_mg, P);
%         fvs_mg = permute(fvs_mg, P);
%         [eps_mg2, ds_mg2] = meshgrid(turb_epss, ds);
%         betas = zeros(size(ds_mg));
%         bs = zeros(length(ds), length(turb_epss));
%         bs_ints_ratio = bs; %zeros(size(ds_mg));
%         %bs_ints_ratio = bs_ints_ratio(:,:,1:length(params.fvs_norm));
% 
%         %Pull numeric densities in this spatial celll
%         z = params.cents_y(iz);
%         cellinds = ((iz-1).*params.Nms + 1):1:(iz.*params.Nms);
%         %Ns = ones(size(cellinds));
%         Ns_cell = ones(size(cellinds));
% 
%         %Pull local fluid properties
%         uL = 0;
% 
%         %Calculate critical breakup diameter
%         d_crit = CalcCritDiameter(iz, params); %m - critical diameter for instability breakage
% 
%         %Create storage matrix
% 
% 
% 
%         %VARY RESOLUTION OF fvs norm based on the number of size groups
%         %below it to reduce computational cost
% 
%         %Iterate through size groups
% 
%         for id = 1:length(ds)
% 
%             %Pull bubble size
%             %mi = params.mms(im); %Current represenative mass
%             %V = (params.nms(im) .* (1 + params.X_mu(iz, im)) .* params.R .* params.T_mu(iz,im))./params.p_func(z);
%             d = ds(id); %(6.*V/pi).^(1/3);
%             %u = params.uzs(cellinds(im));
% 
% 
% 
%             %Iterate through possible epsilons
%             for ie = 1:length(turb_epss)
% 
%                 %Pull epsilon
%                 turb_eps = turb_epss(ie);
%                 params.turb.eps(:) = turb_eps;
% 
%                 %Calculate kolmogorov length scale
%                 lambda_komogorov = ((params.nus.^3)./(turb_eps+1E-8)).^0.25; %m
%                 lambda_min = 31.4 * lambda_komogorov; %Minimum eddy diameter to break a bubble - integrate up to bubble diameter 
% 
%                 %Calculate minimum lambda for breakage
% 
% 
%                 %Calculate breakage rate and BSD
%                 [b_eddy, beta, int_ratio] = BreakageEddyAlt(iz, id, d, Ns_cell, lambda_min, params.break.N_lambdas, params);
% 
%                 %Store results
%                 bs_ints_ratio(id, ie) = int_ratio;
%                 bs(id, ie) = b_eddy; %Breakage rate
%                 betas(id, ie, :) = beta; %Bubble size distribution
% 
% 
%             end
% 
% 
%         end
% 
%         %Create and store interpolation function
%         %betas = permute(betas, P);
%         %F = griddedInterpolant(ds_mg, turb_epss_mg, fvs_mg, betas);
%         params.betas{iz} = @(d, eps, fv) interpn(ds_mg, turb_epss_mg, fvs_mg, betas, d, eps, fv); %function of fv, d, eps
%             %params.bs{iz} = @(d, eps) interp2(ds_mg2, eps_mg2, bs, d, eps);
%         params.beta_ratio{iz} = @(d, eps) interp2(ds_mg2', eps_mg2', bs_ints_ratio', d, eps);
% 
% 
% 
% 
% 
%     end
% 
% 
%     %Test the 
%     if params.debug
% 
%         %Create figure
%         figure();
% 
%         %Reference cases
% 
%         %Test cases
%         eps_ref = 1;
%         ds_ref = [0.0015, 0.002, 0.003, 0.006];
%         BSDs = zeros(length(ds_ref), length(params.fvs_norm_all));
%         for it_debug = 1:length(ds_ref)
% 
%             %Calculated interpolated value
%             BSDs(it_debug,:) = params.betas{iz}(ds_ref(it_debug), eps_ref, params.fvs_norm_all);
% 
%             %Renormalize interpolated value
% 
% 
%             %Plot
%             plot(params.fvs_norm_all, BSDs(it_debug,:), 'LineWidth', 1.5); hold on;
% 
%         end
% 
% 
%         %Plot aesthetics
%         grid on; grid minor; axis square;
% 
% 
%         %Review and close
%         %close;
% 
% 
% 
% 
%     end
% 
%     %If isothermal, convert to standard format for each z-level
%     if params.isothermal
%         for i = 2:params.Nz
%             params.betas{i} = params.betas{1};
%             %params.bs{2:end} = params.bs{1};
%             params.beta_ratio{i} = params.beta_ratio{1};
%         end
%     end
% 
%     %Consider BSDs as a function of d and epsilon for each spatial cell,
%     %which is assumed to have fixed conditions 
% 
%     %Save interpolation function
%     break_out.beta = params.betas;
%     break_out.beta_ratio = params.beta_ratio;
%     outname = sprintf('break_Nd-%d_Nz-%d_Ne-%d_TL-%d', params.break.N_ds, params.Nz, params.break.N_u_spfs, round(params.T_liq));
%     save(outname, 'break_out');
% 
%     %Define outputs
%     params.break.beta_ratio = break_out.beta_ratio;
%     params.break.beta = break_out.beta;
% 
% end


function ub_funcs = BubbleVelocities(params)

    %Calculate range of bubble diameters
    V_min = (params.nms(1) .* params.R .* (params.T_mu_i - 4 * params.T_std_i))./params.ps(1);
    d_min = ((6 .* V_min)./pi).^(1/3);
    V_max = (2 .* params.nms(end) .* params.R .* params.T_liq)./params.p0;
    d_max = ((6 .* V_max)./pi).^(1/3);
    ds = logspace(log10(d_min), log10(d_max), 500);
    
    %For each spatial cell (and its corresponding properties)
    ub_funcs = cell(1, params.Nz);
    for iz = 1:params.Nz
        z = params.zms(iz);
        T = params.T_Lz(z);
        p = params.p_func(z);
        ubs = CalcVelocities(ds, T, params.liquid, p, params.fsolve_opts);
        ub_funcs{iz} = @(d) interp1(ds, ubs, d);
    end




end



function params = LoadTempChars(params)

    if params.heat.active 
        fprintf('-- Generating Characteristics.\n')
        if strcmp(params.chars.type, 'LC')
    
            if params.chars.load
                try
                    if strcmp(params.chars.type, 'LC')
                        filename = sprintf('Chars_%s_Nm=%d_Ti=%d_Tstd=%d_Tl=%d.mat', params.liquid.name, params.Nms,...
                                           round(params.T_mu_i), params.T_std_i, round(params.T_liq)); %NOTE - Must match spec used in TemperatureCharacteristics.m
                        filepath = [params.chars.data_dir, filename];
                        load(filepath);
                        params.T_z.chars = chars;
                    end
                catch
                    warning(sprintf('File (%s) not found, calculating characteristics instead.', filename));
                    params.T_initials = InitializeTemperature(params);
                    params.T_z = TemperatureCharacteristics(params);
                end
    
            else
                params.T_initials = InitializeTemperature(params);
                params.T_z = TemperatureCharacteristics(params);
            end
    
            %params.T_initials = InitializeTemperature(params);
        elseif strcmp(params.chars.type, 'SBM')
            params.T_z = SBM_Characteristics(params); %Load results from SBM studies
        end
    end


end


function params = Loadfvihfuncs(params, disc)

    if params.break.fvih_loadoutput
        try

            %Calculate required filename
            Nz = params.Nz;
            TLtop = params.T_Lz(params.zms(1));
            TLbot = params.T_Lz(params.zms(end));
            filename = sprintf('fv_ih_max_funcs_Nz-%d_TLtop-%d_TLbot-%d.mat', Nz, round(TLtop), round(TLbot));
            load(filename);
            params.fvih_max_funcs = fv_ih_max_funcs;

            %Check if the number of zs is correct
            dataset_N_z = length(params.fvih_max_funcs);
            if ~(dataset_N_z == disc.Nz)
                error('fvih dataset for wrong sized z') %if this occurs, have it generate a new dataset for the current 
            end

        catch
            params.fvih_max_funcs = Calc_fvhi_max_func(params);
        end
    else

        %NOTE - NOTE - Have it store the temperature profile and sigmas to
        %   check if the thing is valid
        % - Change naming scheme to reflect number of z cells
        params.fvih_max_funcs = Calc_fvhi_max_func(params);
    end
end