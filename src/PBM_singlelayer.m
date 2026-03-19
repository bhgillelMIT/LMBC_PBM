function [reactor, disc, sol] = PBM_singlelayer(reactor, disc, sol)



    
    %Only modify settings if there is a single layer
    if sol.single_layer

        %Update user
        fprintf('-- Activating Single Layer Mode.\n')

        %Disable other components
        sol.scheme = 'None';
        sol.orifice_BC_type = 1;
        sol.src_delay = 0;

        %Update discretization parameters
        disc.Nz = 1;

        %Update geometry
        %reactor.H = 0.01;

    end





end