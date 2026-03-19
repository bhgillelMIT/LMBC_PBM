function [Nb_i, Fb_i] = PBM_BCs(reactor, mesh, inlet, disc, zs, rms, Vms, T_orifice, liquid, p_orifice, inputs, fsolve_opts)


    debug = true;

    if disc.Nz == 1

        %Define initial distribution
        u_spf = inputs.reactor.u_spf_orifice;
        L_o = zs(2)-zs(1);
        V_cell_o = reactor.Ac * L_o; %m^3 - Volume of boundary cell(s)
        ub_o = CalcVelocities(2.*rms, T_orifice, liquid, p_orifice, fsolve_opts);
        tres_o = L_o./ub_o; %s - time a bubble spends in cell
        Vres_o = Vms*tres_o; %m3*s

        Fb_i = (1./(inlet.m.std_i .* sqrt(2.*pi))) .* exp(-0.5 .* ((rms - inlet.m.mu_i).^2)./(inlet.m.std_i.^2)); % # of bubbles/m3 - consider integration approach to  calculating
        Nb_i = V_cell_o .* Fb_i; %Initial number of bubbles of each size in the cell


        

        %Calculate
        if inputs.src.solve_alphag
            
            V_dot = sum(Nb_i .* Vms);
            u_spf_o = V_dot/reactor.Ac;
            Nb_scalar = u_spf/u_spf_o;
            Nb_i = Nb_i .* Nb_scalar; Fb_i = Fb_i .* Nb_scalar;
            V_ress = Nb_i .* Vms .* tres_o'; Vres = sum(V_ress);
            alpha_g = Vres/V_cell_o;

        else

            %Scale 
            alpha_g = inputs.src.alphag_manual;
            V_res = V_cell_o * alpha_g; %m3 - Gas volume
            V_res_nom = sum(Nb_i .* Vms);
            V_res_scalar = V_res/V_res_nom;
            Nb_i = Nb_i .* V_res_scalar; Fb_i = Fb_i .* V_res_scalar;

            %Test the volume
            Vs = Nb_i .* Vms;
            V_tot = sum(Vs);
            V_cumsum = cumsum(Vs);
            



        end


        %Normalize to have 1000 bubbles initially
        %Fb_i = 10/sum(Nb_i) .* Fb_i;
        %Nb_i = 10/sum(Nb_i) .* Nb_i;



        %Analytical test cases 
        if strcmp(inputs.src.coalesce_model, 'Scott_1968')
            
            %
            v_bar = 2000 .* disc.V_min;
            v0 = v_bar/2;
            N0 = 1;
            d0 = 0.5 * v_bar;
            vs = Vms;
            vbs = disc.Vbs;
            

            F_func = @(v) (N0/v0) .* (v./v0) .* exp(-v./v0); %(N0/v0) .* (v./v0) .* exp(-v./v0);


            %Calculate integrals
            for i = 1:length(vs)
                vb1 = vbs(i);
                vb2 = vbs(i+1);
                Fb_i(i) = integral(F_func, vb1, vb2);
            end
            %Fb_i = F_func(vs);

            %Fb_i  =  (N0/v0) .* (vs./v0) .* exp(-vs./v0); %(N0/v0) .* exp(-vs./v0);
        
            Fb_i = Fb_i; Nb_i = V_cell_o*Fb_i;

            %Plot the initial distribution
            if debug
                figure();
                loglog(vs, Fb_i);
                ylim([1E-15, 100]); xlim([disc.V_min, disc.V_max])
            end



        elseif strcmp(inputs.src.coalesce_model, 'Hounslow_1988')

            %
            v_bar = 2000 .* disc.V_min;
            v0 = v_bar;
            N0 = 1;
            d0 = 0.5 * v_bar;
            vs = Vms;
            vbs = disc.Vbs;

            F_func = @(v) (N0/v0) .* exp(-v./v0);%(N0/v0) .* (v./v0) .* exp(-v./v0);


            %Calculate integrals
            for i = 1:length(vs)
                vb1 = vbs(i);
                vb2 = vbs(i+1);
                Fb_i(i) = integral(F_func, vb1, vb2);
            end
            %Fb_i = F_func(vs);

            %Fb_i  =  (N0/v0) .* (vs./v0) .* exp(-vs./v0); %(N0/v0) .* exp(-vs./v0);
        
            Fb_i = Fb_i; Nb_i = V_cell_o*Fb_i;

            %Plot the initial distribution
            if debug
                figure();
                loglog(vs, Fb_i);
                ylim([1E-15, 100]); xlim([disc.V_min, disc.V_max])
            end

        elseif strcmp(inputs.src.breakage_model, 'Uniform_Binary')
            
            Fb_i = zeros(size(Fb_i));
            [nearestval, nearestind] = min(abs(Vms - 1));
            Fb_i(nearestind) = 1/(disc.Vbs(nearestind+1) - disc.Vbs(nearestind));



        end

        
        

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