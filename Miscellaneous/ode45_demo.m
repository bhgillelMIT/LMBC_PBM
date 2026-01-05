function ode45_demo()

%% Setup

    %Plot settings
    lw = 2;
    fs = 18;

%% Model
    
    %Inputs
    ti = 0; %s - initial time
    tf = 10; %s - stop time
    dt = 0.01; %s - timestep
    zi = 100; %m - initial height
    ui = 0; %m/s - initial velocity
    params.d = 0.001; %m - droplet diameter
    params.g = 9.81; %m/s2 - gravitational acceleration
    params.rho_l = 6500; %kg/m3 - density of the droplet
    params.rho_g = 1; %kg/m3 - density of gas
    params.Cd = 0.5; %unitless - drag coefficient of droplet
    params.ug = 0; %m/s - gas velocity; positive = upwards
    params.dt_max = dt;

    %Calculate droplet initial properties
    params.V = (4/3) .* pi .* (params.d/2).^3; %m3 - droplet volume
    params.Ac = pi .* (params.d/2).^2;
    params.As = 4 * params.Ac;
    params.m = params.V .* params.rho_l; %kg - droplet mass

    %Define ode function
    ode_func = @(t,y) falling_droplet(t,y,params);

    %Define initial condition
    y0(1) = zi;
    y0(2) = ui;

    %Define ode options
    odeopts = odeset('InitialStep', 1e-3, 'MaxStep', params.dt_max, 'RelTol', 1e-6, 'AbsTol', 1e-6, 'Stats','on'); 

    %Solve system
    tspan = ti:dt:tf;
    [T,Y] = ode45(ode_func,tspan, y0, odeopts);

    %Plot results
    figure();
    subplot(1,2,1);
    plot(T, Y(:,1), 'k-', 'LineWidth', lw);
    xlabel('Time (s)'); ylabel('Z-position (m)');
    title('Z-position');
    grid on; grid minor; axis square;
    set(gca, 'FontSize', fs);

    subplot(1,2,2);
    plot(T, Y(:,2), 'b-', 'LineWidth', lw);
    xlabel('Time (s)'); ylabel('Z-velocity (m/s)');
    title('Z-velocity');
    grid on; grid minor; axis square;
    set(gca, 'FontSize', fs);


end


function dydt = falling_droplet(t,y,params)

    %Isolate variables
    z = y(1); %z is relative to a fixed reference frame
    u = y(2); 

    %Calculate terminal velocity - for reference 
    u_term = sqrt((2.*params.m.*params.g)./(params.Cd .* params.rho_g .* params.Ac));


    %Calculate force balance on droplet
    Fg = params.m .* params.g; %N - gravitational force 
    Fd = 0.5 .* params.Cd .* params.rho_g .* (u+params.ug).^2 .* params.Ac; %N - drag force

    %Define physical terminating condition
    if z > 0 
        dzdt = y(2)+params.ug; %Time derivative of height is just the current velocity
        dudt = (Fd - Fg)/params.m;
    else %Droplet no longer falling
        dzdt = 0;
        dudt = 0;
    end

    %Define output - an array of first order derivatives. For higher order
    %derivatives, express them as a collection of first order derivatives.
    %For instance: d2zdt2 = d/dt(dz/dt)
    dydt = [dzdt; dudt]; %Order must match the order above - Must be vertical vector

    


end