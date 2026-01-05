function [Nb_i, Fb_i] = PBM_BCs(reactor, mesh, inlet, disc, zs, rms, Vms, T_orifice, liquid, p_orifice, fsolve_opts)

    if disc.Nz == 1

        %Define initial distribution
        L_o = zs(2)-zs(1);
        V_cell_o = reactor.Ac * L_o; %m^3 - Volume of boundary cell(s)
        ub_o = CalcVelocities(2.*rms, T_orifice, liquid, p_orifice, fsolve_opts);
        tres_o = L_o./ub_o; %s - time a bubble spends in cell
        Fb_i = (1./(inlet.m.std_i .* sqrt(2.*pi))) .* exp(-0.5 .* ((rms - inlet.m.mu_i).^2)./(inlet.m.std_i.^2)); % # of bubbles/m3 - consider integration approach to  calculating
        Nb_i = V_cell_o .* Fb_i; %Initial number of bubbles of each size in the cell
        
        %Normalize to have 1000 bubbles initially
        Fb_i = 10/sum(Nb_i) .* Fb_i;
        Nb_i = 10/sum(Nb_i) .* Nb_i;
        

    else
    
        %Define initial distribution
        L_o = zs(2)-zs(1);
        V_cell_o = reactor.Ac * L_o; %m^3 - Volume of boundary cell(s)
        ub_o = CalcVelocities(2.*rms, T_orifice, liquid, p_orifice, fsolve_opts);
        tres_o = L_o./ub_o; %s - time a bubble spends in cell
        Fb_i = (1./(inlet.m.std_i .* sqrt(2.*pi))) .* exp(-0.5 .* ((rms - inlet.m.mu_i).^2)./(inlet.m.std_i.^2)); % # of bubbles/m3 - consider integration approach to  calculating
        Nb_i = V_cell_o .* Fb_i; %Initial number of bubbles of each size in the cell
        Nb_i_leaving = Nb_i .* 1; % 1/s - number of bubbles leaving (and entering) per second
    
        Vdot_i = (Nb_i .* Vms)./tres_o'; %m3/s 
        Vdot_i_sum = sum(Vdot_i);
        N_bubbles_s = reactor.V_dot./Vdot_i_sum;
    
       % Vdot_check = Vms(2) * N_bubbles_s;
    
        %Calculate flux
        gauss_func = @(r) (1./(inlet.m.std_i .* sqrt(2.*pi))) .* exp(-0.5 .* ((r - inlet.m.mu_i).^2)./(inlet.m.std_i.^2));
        Nb_i = gauss_func(rms); %1/s - number of bubbles crossing border
        V_dot_i = Nb_i .* Vms; %m3/s
        V_dot = sum(V_dot_i);
        V_dot_scalar = reactor.V_dot./V_dot;
        Nb_i = Nb_i .* V_dot_scalar; %Initial
    
        %Check
        V_dot_check = sum(Nb_i .*  Vms);
        u_spf_check = V_dot_check/reactor.Ac;
    
         %Vb_i = Vms .* Fb_i .* V_cell_o;
        %Vb_i_sum = sum(Vb_i);
        %N_bubbles_s = reactor.V_dot/Vb_i_sum; %Number of bubbles being generated - multiply by distribution
        Fb_i = N_bubbles_s * Fb_i;
        Fb_i(Fb_i < 0.1) = 0;
        Nb_if = Fb_i * V_cell_o;
        Vb_if = Nb_if .* Vms;
        eps_gas_i = sum(Vb_if)/V_cell_o;
        
    
    
        %Estimate the gas fraction in the first cell
        V_gas = Nb_i .* (mesh.hs(1)./ub_o') .* Vms;
        V_gas = sum(V_gas);
        alpha_gas = V_gas./V_cell_o;
    end




end