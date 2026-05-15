% This function performs a convergence study on the Population Balance Model 
% by considering several different mesh resolutions

function PBM_convergence()

%% Create water structure

    water.dyn_visc = @(T) 0.0010016 .* ones(size(T)); %Pa*s - at 20 C
    water.mol_mass = 0.0180153; %kg/mol
    water.surf_tension =  @(T) 0.0728 .*ones(size(T)); 
    water.density = @(T) 998.21 .* ones(size(T)); 
    water.Cp = @(T) 4157 .* ones(size(T)); 
    water.therm_cond = @(T) 0.598 .* ones(size(T)); 
    water.name = 'water';

%% Convergence Study

    %Specify range of mesh sizes to test
    N_min = 5
    N_max = 200;
    N_sizes = 30;

    %Specify test case conditions
    u_spf = 0.05;
    epss = 10 .* u_spf;
    alpha_gs = 0.154;


    %Define mesh sizes
    N_pts = round(logspace(log10(N_min), log10(N_max), N_sizes));

    
    inputs.sol.src_delay = 0;
    inputs.sol.break_file = 'break_water_Nd-20_Nz-1_Ne-20_TL-293.mat';
    inputs.sol.postprocess = false; %Disabled to prevent figures 

    %Initialize results structure
    results = struct('N', cell(1, N_sizes), 'y_final', cell(1, N_sizes));

    %Iterate through cases
    for in = 1:N_sizes



        %Extract current mesh size
        N = N_pts(in);
        fprintf('Running case with N = %d\n', N);

        %Initialize inputs 
        inputs = PBM_inputs();
        inputs.disc.Nms = N;
        inputs.disc.mesh_hybrid_cells = round(inputs.disc.Nms .* 0.5);  
        inputs.reactor.liquid = water;
        inputs.reactor.u_spf_orifice = u_spf; %m/s
        inputs.reactor.T_gas_i_mu = 20 + 273.15;
        inputs.reactor.T = 20 + 273.15;
        inputs.reactor.T_orifice = inputs.reactor.T;
        inputs.reactor.T_min_i = 0 + 273.15;
        inputs.src.eps_manual = epss(ic);
        inputs.src.solve_eps = false;
        inputs.src.alphag_manual = alpha_gs(ic);
        inputs.src.solve_alphag = false;
        inputs.sol.heat.active = false;
        inputs.sol.react.active = false;
        inputs.sol.solve_ub = false;
        inputs.sol.single_layer = true;
        inputs.sol.ub_manual = inputs.reactor.u_spf_orifice/inputs.src.alphag_manual; 
        inputs.reactor.H = 0.3;
        inputs.src.breakage_active = true;
        inputs.src.c_wake = Kws(ik);
        inputs.sim.t_end = 1;

    
        
        %Run PBM with current mesh size
        output = PBM_v3(N);
        
        %Extract solution at final time
        y_final = output.Y(end, :);
        
        %Store results for convergence analysis
        results(in).N = N;
        results(in).y_final = y_final;
    end




end