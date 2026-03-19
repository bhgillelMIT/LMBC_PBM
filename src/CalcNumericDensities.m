function [Ns_z, Ns_m, Ns_T, Ns_fracs] = CalcNumericDensities(y, params)

    %Allocate output vector
    Ns_m = params.Ns_m_zero;
    Ns_z = params.Ns_z_zero;
    Ns_T = params.Ns_T_zero;

    %Iterate through each spatial cell
    ia = 1;
    iT = 1;
    inds = [];
    for iz = 1:params.Nz
        for im = 1:params.Nms
            
            %Iterate through temperature bins
            minds = find(params.zinds == iz & params.xinds == im);
            Tinds = unique(params.Tinds(minds));
            for it = 1:length(Tinds)

                %Iterate through conversion bins
                Xinds = unique(params.Xinds(minds));
                for ix = 1:length(Xinds)

                    %Find overall index
                    ind = find(params.zinds == iz & params.xinds == im & params.Tinds == Tinds(it) & params.Xinds == Xinds(ix));

                    %Add the numeric density
                    Ns_z(iz) = Ns_z(iz) + y(ind);
                    Ns_m(ia) = Ns_m(ia) + y(ind);
                    Ns_T(iT) = Ns_T(iT) + y(ind);

                    %inds(end+1) = ind;

                end

                %Update iteration counter 
                iT = iT + 1;

            end

            %Update iteration counter
            ia = ia + 1;

        end
    end

    %Calculate the mass fraction of each cell within each representative
    %size
    it = 1;
    Ns_fracs = zeros(size(y));
    for iz = 1:params.Nz
        for im = 1:params.Nms

            %Pull total numeric density of this representative size
            Nm = Ns_m(it);

            %Identify sub-indices
            subinds = find(params.zinds == iz & params.xinds == im);

            %Calculate mass fraction of each
            for sind = subinds
                if Nm > 0
                    Ns_fracs(sind) = y(sind)/Nm;
                else
                    Ns_fracs(sind) = 0;
                end
            end

            %Update iteration counter
            it = it + 1;


        end
    end


end