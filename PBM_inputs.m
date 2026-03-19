function inputs = PBM_inputs()
% This function outputs the default structure for PBM inputs. It should be
% edited to change default settings for the PBM.

    %Material Properties
    addpath('Material Properties/');
    load tin.mat
    load hydrogen.mat
    load carbon.mat
    load methane.mat
    load argon.mat
    load water.mat

    %Add name to tin
    tin.name = 'tin';
    water.name = 'water';

    %Main settings
    sim.t_end = 10;
    reactor.u_spf_orifice = 0.0037; %m/s
    inlet.m.mu_i = 0.006; %0.003;
    inlet.m.std_i = 0.00075; %%0.00075;

    %Animation settings
    animate.active = true;
    animate.save = true;
    animate.fname = 'PBM_gif_coalesce_v3.gif';
    animate.gifdelay = 0.25;
    animate.time_dilation = 1;
    inputs.fsolve_opts = optimoptions('fsolve', 'Display', 'off');

    %Simulation
    sol.src_delay = 0;
    
    

    %Geometry
    reactor.D = 0.15;
    reactor.R = reactor.D./2;
    reactor.H = 0.3;
    reactor.Ac = pi .* reactor.R^2;

    %Operating parameters
    reactor.T = 293.15; %1300 + 273.15;
    reactor.T_orifice = reactor.T;
    reactor.T_min = 273.15; %800 + 273.15;
    reactor.X_Ar = 0.9;
    
    reactor.p_surf = 101325;
    reactor.M_gas = 0.016;
    reactor.liquid = water;
    reactor.T_gas_i_mu = 293.15; %1000 + 273.15;
    reactor.T_gas_i_std = 2; %15;
    reactor.T_min_i = reactor.T_min; %800 + 273.15; %Minimum temperature to consider for heat transfer
    

    % %Calculate reactor floew rate
    % reactor.V_dot = reactor.u_spf_orifice .* reactor.Ac; %m3/s
    % reactor.V_dot_Lmin = 1000 .* reactor.V_dot .* 60;


    %Discretization
    mmesh.type = 'Geometric'; %Type of mesh to use. Options: 'Linear', 'Geometric', 'Hybrid'
    mmesh.hybrid_cells_frac = 0.5;
    mmesh.input = 'Radius'; %Options: 'Radius', 'Volume'
    disc.Nr = 1;
    disc.Nz = 10; 
    disc.Nms = 50;
    disc.Nrs = disc.Nms; % # of discrete volumes considered
    disc.mesh_hybrid_cells = round(disc.Nms .* mmesh.hybrid_cells_frac);  
    disc.r_min = 0.5 .* 0.001; %m - smallest size of bubble to consider
    disc.r_max = 15 .* 0.001;
    disc.dt = 0.01; %initial step size, solver is adaptive
    disc.hybrid_cells_frac = mmesh.hybrid_cells_frac;
    %disc.mesh_hybrid_frac = mmesh.hybrid_frac;

    %Initial condition
    initial.u = 0;
    initial.T = 573.15;

    %Source term controls
    src.solve_eps = true;
    src.eps_manual = 0.5;
    src.solve_alphag = true;
    src.alphag_manual = 0.154;
    src.coalesce_model = 'Wang_2005'; %Options: 'Wang_2005', PrinceBlanch_1990', 'Scott_1968', 'Hounslow_1988' 
    src.breakage_model = 'Wang_2005'; %Options: 'Wang_2005', 'Uniform_Binary', 'Luo_Svendson_1996'
    src.coalesce_constant_rate = 1;
    src.breakage_constant_rate = 1;
    src.coalesce_active = true;
    src.breakage_active = true;
    src.radial_active = true;
    src.radial_alpha_g_func = @(x) 
    src.radial_ul_fnc
   

    %Soler controls
    sol.type = 'direct';
    sol.scheme = 'Upwind';
    sol.heat.active = false; %boolean - Determines if heat transfer is considered in the simulation
    sol.react.active = false; %boolean - Determines if the reaction is considered in the simulation
    sol.solve_ub = true; %boolean - Determines if the bubble velocity is solved each iteration, or imposed initially - used only for test cases
    sol.ub_manual = 0.3; %m/s - Velocity of bubbles; Must either be a single value (in which case all bubbles travel at the same velocity), or a vector with a length matching the number of representative masses (params.Nms)
    sol.solve_detail_its = 1; %iterations - number of iterations between calculations of certain quantities in the simulation that have significant costs  associated with them. Use the most recent value 
    sol.disable_advection = false;
    sol.disable_advection_t = 10;
    sol.single_layer = true ;
    sol.orifice_BC_type = 2;
    sol.break_file = 'break_tin_Nd-25_Nz-1_Ne-25_TL-1573.mat';

    %Boundary conditions
    

    %Define input structure
    inputs.animate = animate;
    inputs.reactor = reactor;
    inputs.disc = disc;
    inputs.initial = initial;
    inputs.inlet = inlet;
    inputs.mmesh = mmesh;
    inputs.sim = sim;
    inputs.src = src;
    inputs.sol = sol;

end