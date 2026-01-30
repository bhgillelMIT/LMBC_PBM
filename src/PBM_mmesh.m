function [disc, mmesh] = PBM_mmesh(disc, mmesh, reactor, p_orifice, rho_gas_orifice, debug)

    R = 8.3145;

    %Calculate mass - case with mass as size dimension
    if strcmp(mmesh.input, 'Radius')

        disc.V_min = (4/3) .* pi .* disc.r_min^3;
        disc.V_max = (4/3) .* pi .* disc.r_max^3;
    elseif strcmp(mmesh.input, 'Volume')
        disc.V_min = disc.r_min;
        disc.V_max = disc.r_max;
    else
        error('Invalid input type for mmesh.input; Please specify either "Radius" or "Volume".');
    end
    disc.m_min = disc.V_min * rho_gas_orifice;
    disc.m_max = disc.V_max * rho_gas_orifice;
    switch mmesh.type
        case 'Geometric'
            disc.m_rat = (disc.m_max/disc.m_min).^(1./disc.Nms);
            disc.mbs = disc.m_min .* disc.m_rat.^(0:1:disc.Nms);
        case 'Linear'
            disc.mbs = linspace(disc.m_min, disc.m_max, disc.Nms+1);
        case 'Hybrid'

            %Calculate cutoff 
            disc.N_linear = disc.Nms - disc.mesh_hybrid_cells;
            disc.N_logari = disc.mesh_hybrid_cells;

            X_low = 0;
            X_high = 1;
            X_mid = 0.5;
            err = 1;
            while err > 1E-12

                %Calculate m_mid
                disc.m_mid = (disc.m_max - disc.m_min) * X_mid + disc.m_min;

                %Solve logarithmic section
                disc.m_rat = (disc.m_mid/disc.m_min).^(1./(disc.mesh_hybrid_cells+1));
                disc.mbs_small = disc.m_min .* disc.m_rat.^(1:disc.mesh_hybrid_cells+2);

                disc.dms_small = diff(disc.mbs_small);

                %Solve linear section
                disc.mbs_large = linspace(disc.mbs_small(end-1),disc.m_max, disc.N_linear + 1);
                disc.mbs_large = disc.mbs_large(2:end);
                disc.dms_large = diff(disc.mbs_large);

                %Compile results
                disc.mbs = [disc.mbs_small(1:end-1), disc.mbs_large];

                %Define error and determine how to change X_mid
                err = disc.dms_large(1) - disc.dms_small(end);
                if err > 0
                    X_low = X_mid;
                    X_mid = 0.5 * (X_high - X_mid) + X_mid;   
                else
                    X_high = X_mid;
                    X_mid = 0.5 * (X_mid - X_low) + X_low;
                end
                err = abs(err);


            end

            %Debug plot
            if debug
                figure(); 
                semilogx(disc.mbs, zeros(size(disc.mbs)), 'r.')
                xlabel('Mass (kg)');
            end




            % disc.mesh_hybrid_frac = 1;
            % 
            % 
            % %Fill bins
            % 
            % disc.m_mid = disc.m_min + (disc.m_max - disc.m_min) .* disc.mesh_hybrid_frac;
            % disc.m_rat = (disc.m_mid/disc.m_min).^(1./(disc.mesh_hybrid_cells));
            % disc.mbs_small = disc.m_min .* disc.m_rat.^(1:disc.mesh_hybrid_cells);
            % disc.mbs_large = linspace(disc.mbs_small(end),disc.m_max, disc.Nms-disc.mesh_hybrid_cells+2);
            % disc.mbs_large = disc.mbs_large(2:end);
            % 
            % % disc.mbs_small = %linspace(disc.m_min, disc.m_mid, disc.mesh_hybrid_cells);
            % % disc.m_rat = (disc.m_max/disc.mbs_small(end)).^(1./(disc.Nms-disc.mesh_hybrid_cells));
            % % disc.mbs_large = disc.mbs_small(end) .* disc.m_rat.^(1:(disc.Nms-disc.mesh_hybrid_cells+1));
            % disc.mbs = [disc.mbs_small, disc.mbs_large];

    end
    mbs = disc.mbs;
    mms = (mbs(1:end-1) + mbs(2:end))./2; %(mbs(1:end-1) + diff(mbs))./2;
    rbs = ((3*mbs*R*reactor.T_orifice)./(4*pi*reactor.M_gas*p_orifice)).^(1/3);
    rms = ((3*mms*R*reactor.T_orifice)./(4*pi*reactor.M_gas*p_orifice)).^(1/3);
    Vms = (mms .* R .* reactor.T_gas_i_mu)./(p_orifice .* reactor.M_gas); %m3
    Vbs = (mbs .* R .* reactor.T_gas_i_mu)./(p_orifice .* reactor.M_gas);
    nms = (Vms * p_orifice)./(R .* reactor.T_gas_i_mu); %mols
    nbs = (Vbs * p_orifice)./(R .* reactor.T_gas_i_mu);

    %Calculate radii - case with radius as size dimension
    disc.rb_brackets = linspace(disc.r_min, disc.r_max, disc.Nrs+1); %logspace(log10(r_min), log10(r_max), Nrs+1);
    disc.rb_mids = disc.rb_brackets(1:end-1) + diff(disc.rb_brackets)./2;
    disc.Vbs = Vbs;

    %Define mmesh structure
    mmesh.mbs = mbs;
    mmesh.mms = mms;
    mmesh.rbs = rbs;
    mmesh.rms = rms;
    mmesh.Vbs = Vbs;
    mmesh.Vms = Vms;
    mmesh.nms = nms;
    mmesh.nbs = nbs;

end