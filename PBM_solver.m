function  [T,Y] = PBM_solver(tspan, y0, params, options, varargin)

%% Setup 


    %Make parameters global
    global params
    
    solver_name = 'PBM_solver';

    %Define ODEs
    %ode = @(t,y) PBM_ode(t,y,params);

%% Process Inputs

    % Handle fewer inputs
    if nargin < 4
        options = [];
        if nargin < 3
            y0 = [];
            if nargin < 2
                tspan = [];
            end
        end
    end

    % % Stats
    % nsteps   = 0;
    % nfailed  = 0;
    % nfevals  = 0;
    % npds     = 0;
    % ndecomps = 0;
    % nsolves  = 0;

    % % Repackage ODE
    % [ode, odeIsFuncHandle, odeTreatAsMFile] = packageAsFuncHandle(ode);
    % 
    % % Output
    % output_sol = odeIsFuncHandle && nargout == 1;   % sol = odeXX(...)
    % output_ty  = ~output_sol && nargout > 0;        % [t,y,...] = odeXX(...)
    % % There might be no output requested...
    % 
    % sol = []; kvec = []; dif3d = [];
    % if output_sol
    %     sol.solver = solver_name;
    %     sol.extdata.odefun = ode;
    %     sol.extdata.options = options;
    %     sol.extdata.varargin = varargin;
    % end
    % 
    % % Handle solver arguments
    % [neq, tspan, ntspan, next, t0, tfinal, tdir, y0, f0, odeArgs, ...
    %     options, threshold, rtol, normcontrol, normy, userhmin, hmax, htry, htspan] = ...
    %     odearguments(odeIsFuncHandle, odeTreatAsMFile, solver_name, ode, tspan, y0, options, varargin);
    % nfevals = nfevals + 1;



    %Settings
    ti = tspan(1);
    tf = tspan(end);
    dt_i = 0.025; %s

    %Specify options for ode solvers
    ode45opts = odeset('InitialStep', 1e-3, 'MaxStep', params.dt_max, 'RelTol', 1e-4, 'AbsTol', 1e-6, 'Stats','on', 'OutputFcn', @PBM_output);
    ode15opts = odeset('InitialStep', 1e-3, 'MaxStep', params.dt_max, 'RelTol', 1e-4, 'AbsTol', 1e-6, 'Stats','on', 'OutputFcn', @PBM_output);




    switch params.sol.type
        case 'segregated'
    
            %Print update
            fprintf('PBM - Using Segregated Solver w/ Strang Splitting.\n\n');

            %Isolate settings
            sep_layer = params.sol.sep_layer;

            % 
            % 
            % %Calculate Jacobian for the ode system
            % [Jconstant,Jac,Jargs,Joptions] = ...
            %     odejacobian(odeIsFuncHandle,ode,t0,y0,options,varargin);
            % Janalytic = isempty(Joptions);
        
            
            
        
        
            %Allocate storage
            N_steps = ceil(tf/dt_i);
            Yadv = zeros(2*N_steps+1, length(y0));
            Ysrc = zeros(N_steps, length(y0));
            Yfin = zeros(N_steps, length(y0));
            Tfin = zeros(N_steps, 1);
        
            %Set initial values
            Yadv(1,:) = y0'; Yadv(2,:) = y0';
            Ysrc(1,:) = y0';
            Yfin(1,:) = y0';
            Tfin(1) = ti;
        
        
            it = 1;
            t = ti;
            dt = dt_i;
            dth = 0.5*dt;
            while (t + dt) < tf
        
                %Update user
                fprintf('PBM (t = %0.6f s; it = %d', t, it)
        
                %Handle non-uniform final timestep
                if (t + dt) > tf
                    dti = tf - t;
                    dthi = dti/2;
                else
                    dti = dt;
                    dthi = dth;
                end
        
                %Half Advection
                tspan = [t, t+dth];
                params.sol.sep_layer = false; %Turn off for advection
                [T,Y1] = ode45(@(t,y) PBM_Advection(t,y,params), tspan, y0, ode45opts);
                y0 = Y1(end,:);
                fprintf('- Advection Step 1 Done;\n');
        
        
                %Full Coalescence and  Breakage -  solve each layer independently
                tspan = [t, t+dt];

                if sep_layer
                    params.sol.sep_layer = sep_layer;
                    for iz = 1:params.Nz

                        params.iz = iz;
                        [T,Y2] = ode45(@(t,y) PBM_Source(t,y,params), tspan, y0, ode45opts);
                        y0 = Y2(end,:); %

                        %Debug
                        if params.sol.debug
                            dY = y0 - Y1(end,:);           
                        end

                    end
                else
                    [T,Y2] = ode45(@(t,y) PBM_Source(t,y,params), tspan, y0, ode45opts);
                end
                fprintf('-- Source Step Done;\n');
                
                
                %Half Advection
                params.sol.sep_layer = false; %Turn off for advection
                tspan = [t+dth, t+dt]; y0 = Y2(end,:);
                [T,Y3] = ode45(@(t,y) PBM_Advection(t,y,params), tspan, y0, ode45opts);
                fprintf('-- Advection Step 2 Done;\n');


                %Define y0 for next step
                y0 = Y3(end,:);
        
                %Log timestep
                Yadv(2*it, :) = Y1(end,:);
                Yadv(2*it+1, :) = Y3(end,:);
               
                
        
        
                %Update time and iteration counter
                t  = t + dt;
                it = it + 1;
        
                %Log data
                Ysrc(it,:) = Y2(end,:);
                Yfin(it, :) = Y3(end,:);
                Tfin(it) = t;
            
            end
        
            %Define outputs
            T = Tfin;
            Y = Yfin;

        otherwise % use solver directly 

            %Print update
            fprintf('PBM - Using Direct Solver (ode15s).\n\n');

            %Initialize with ode45 to improve accuracy of Jacobian
            tspan1 = [tspan(1), 0.05]; %min([0.1, 0.025*tspan(end)])];
            [T1,Y1] = ode45(@(t,y) PBM_ode(t,y,params), tspan1, y0, ode45opts);
            ms1 = params.m_total;

            %Full run
            y0 = Y1(end,:);
            tspan2 = [tspan1(end), tspan(end)];
            [T,Y] = ode15s(@(t,y) PBM_ode(t,y,params), tspan2, y0, ode15opts);
            ms2 = params.m_total;

            %Consolidate outputs
            T = [T1; T]; Y = [Y1; Y];
            ms_total = [ms1(ms1 > 0), ms2(ms2 > 0)];
            params.m_total = ms_total;
            x = 1;



    end




end