%A function to identify cells that match THis should be run in PBM_v3 before running the simulation.

function params = IdentifyAdjacentInds(params, N_volumes)

    %Allocate storage vectors
    params.inds_below = cell(1, N_volumes);
    params.inds_above = cell(1, N_volumes);

    %Iterate through volumes
    for i = 1:N_volumes

        %Pull indexes
        xind = params.xinds(i);
        zind = params.zinds(i);
        Tind = params.Tinds(i);

        if params.heat.active

            %Calculate counts
            NT = unique(params.Tinds(params.xinds == xind & params.zinds == zind));
            NX = 1; %TO BE UPDATED
    
            %Pull temperature characteristics
            merged = params.T_z.chars.merged{xind};

            %Flip merged bins to simplify indexing
            bins_in = flipud(merged.bins_in);
    
            %Identify index below
            if zind == 1
                params.inds_below{i} = [-1];
            else
                
                %Identify bins feeding into this bin
                row_inds = find(bins_in(:,zind) == Tind);
                Tinds_below = bins_in(row_inds, zind-1);

                %Convert into overall index (account for space/mass)
                zind_below = zind - 1;
                inds_below = zeros(size(Tinds_below));
                for ii = 1:length(Tinds_below)
                    inds_below(ii) = find(params.zinds == zind_below & params.xinds == xind...
                    & params.Tinds == Tinds_below(ii));
                end

                inds_below_unique = unique(inds_below);
                % if length(inds_below_unique) ~= length(inds_below)
                %     x = 1;
                % end

                %Store result
                params.inds_below{i} = inds_below_unique;
  
            end
    
            %Identify index above 
            if zind == params.Nz
                params.inds_above{i} = -1;
            else

                %Identify temperature index
                row_inds = find(bins_in(:,zind) == Tind);
                Tinds_above = unique(bins_in(row_inds,zind+1));

                %Convert into overall index (account for space/mass)
                zind_above = zind + 1;
                inds_above = zeros(size(Tinds_above));
                for ii = 1:length(Tinds_above)
                    inds_above(ii) = find(params.zinds == zind_above & params.xinds == xind...
                        & params.Tinds == Tinds_above(ii));
                end
                
                %Store result
                params.inds_above{i} = inds_above;
            end
    
            %params.inds_above
            %params.inds_below


        else

            %Below index
            if zind == 1
                params.inds_below{i} = -1;
            else 
                params.inds_below{i} = i - params.Nms;
            end

            %Above index
            if zind == params.Nz
                params.inds_above{i} = -1;
            else
                params.inds_above{i} = i + params.Nms;
            end

        end


    end
    




end