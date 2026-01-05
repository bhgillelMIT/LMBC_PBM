function dydt = BubbleRise_odefun(t,y, ps)

    %Refernece equations
    % u(0) = 0;
    % z(0) = 0;
    % V(0) = nRT/p = mRT/pM
    % p = p0 + rho*g*(h-z)
    % V = (4/3) * pi * r^3 -> r = (3V/(4*pi))^(1/3);
    % As = 4*pi*r^2 = 4 * pi * 
    % dVdt = nRt * d/dt(1/p) = 
    % dudt = (F_B - F_D)/m = ((rho_L * V * g - m)*g -
    %                           0.5*rho_L*Cd*AAs*u^2)/m
    % dzdt = u
    % dpdt = -rho*g*dzdt
    % Vector = [z, u, p]
    % 

    %Extract parameters
    R = ps.R;
    Cd = ps.Cd;
    %As = ps.As;
    r0 = ps.r0;
    m = ps.m;
    p0 = ps.p0;
    h = ps.h;
    rho_L = ps.rho_L;
    g = ps.g;
    M = ps.M;
    T = ps.T;

%% Algebraic Equations

    
    
    V = (m*R*T)./(M*y(3));
    r = (3*V/(4*pi))^(1/3);
    Ac = pi * r.^2;
    m_vm = (11/16) .* rho_L .* V;

%% Differential Equations

    dydt = zeros(3,1);
    dydt(1) = y(2);
    dydt(2) = ((rho_L * V - m)*g - 0.5*rho_L*Cd*Ac*y(2)^2)/m_vm;
    dydt(3) = -rho_L*g*dydt(1);

    %Track other values


end