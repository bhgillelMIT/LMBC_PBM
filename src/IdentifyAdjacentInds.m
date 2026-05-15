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
        Xind = params.Xinds(i);

        


        %Handle cases
        if params.heat.active & params.react.active

            itrails = params.T_z.chars.itrails;

            %Determine if this is a trailing bin
            trail_inds = itrails{xind};
            trail = any(trail_inds == trail_inds);

            %Calculate counts
            NT = unique(params.Tinds(params.xinds == xind & params.zinds == zind));
            NX = params.chars.N_Xis; %TO BE UPDATED
    
            %Determine if this is a trailing ind 


            %Iterate through conversion

            %Pull temperature characteristics
            merged = params.T_z.chars.merged{xind}{Xind};

            %Flip merged bins to simplify indexing
            bins_in = flipud(merged.bins_in);

            %Pull the conversion characteristics
            Xs_bs = params.T_z.chars.Xs_bs{xind}{Tind}{Xind}; %Conversions of this characteristic at bin boundaries (top/bottom)
            X_in = Xs_bs(zind); %Conversion entering
            X_out = Xs_bs(zind+1); %Conversion exiting 
    
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

                    %Identify the conversion bins above this one
                    XTinds_below = find(params.zinds == zind_below & params.xinds == xind);
                    Xinds_below = find(params.zinds == zind_below & params.xinds == xind...
                        & params.Tinds == Tinds_below(ii));


                    %Pull boundary values
                    Tc_below = zeros(size(Xinds_below));
                    Tb_below = Tc_below;
                    Xc_below = Tc_below;
                    Xb_below = Tc_below;
                    Xb_final = Tc_below;
                    for ia = 1:length(Xinds_below)

                        %Pull indices
                        XTind = XTinds_below(ia);
                        ix = params.Xinds(XTind);
                        it = params.Tinds(XTind);

                        %Determine the temperatures above
                        Tcs_below = params.T_z.chars.merged{xind}{ix}.Tcz(:,zind_below);
                        active_inds = Tcs_below > 0;
                        Tcs_below = Tcs_below(active_inds);
                        Tc_below(ia) = Tcs_below(it);

                        %Determine the conversions above
                        Xcs_below = params.T_z.chars.merged{xind}{ix}.Xcz(:,zind_below);
                        Xcs_below = Xcs_below(active_inds);
                        Xc_below(ia) = Xcs_below(it);
                        Xbs_below = params.T_z.chars.merged{xind}{ix}.Xbz(:,zind_below);
                        Xbs_below = Xbs_below(active_inds);
                        Xb_below(ia) = Xbs_below(it);
                        Xbs_final = params.T_z.chars.merged{xind}{ix}.Xbz(:,end);
                        Xbs_final = Xbs_final(active_inds);
                        Xb_final(ia) = Xbs_final(it);

                    end

                    %Determine bin above 
                    if length(Xinds_below) == params.chars.N_Xis %Simple case - direct mapping
                        inds_below(ii) = find(params.zinds == zind_below & params.xinds == xind...
                        & params.Tinds == Tinds_below(ii) & params.Xinds == Xind);
                    else
                        Xinds_X_below = params.Xinds(Xinds_below);
                        inds_X_equal = Xb_below == X_out;
                        [~,inds_X_closest] = min(abs(Xb_below - X_out));
                        if any(Xinds_X_below == Xind & inds_X_equal) 
                            inds_same = Xinds_below(Xinds_X_below == Xind);
                            inds_below(ii) = inds_same;
                        elseif all(Xb_above < X_out)
                            inds_below(ii) = Xinds_below(end);
                        elseif all(Xb_above > X_out)
                            inds_below(ii) = Xinds_below(1);
                        else
                            inds_below(ii) = Xinds_below(inds_X_closest);
                        end
                       
                    end



                    % 
                    % 
                    % %Pull Conversion values at boundaries
                    % Xs_bs_below = zeros(params.Nz+1, params.chars.N_Xis);
                    % for ix = 1:params.chars.N_Xis
                    %     Xs_bs_below(:,ix) = params.T_z.chars.Xs_bs{xind}{Tinds_below(ii)}{ix};
                    % end
                    % [Xs_bs_f, Xs_bs_order] = sort(Xs_bs_below(end,:));
                    % Xs_bs_below_orig = Xs_bs_below;
                    % Xs_bs_below = Xs_bs_below(:, Xs_bs_order);
                    % Xs_bs_below_layer = Xs_bs_below(zind,:);
                    % 
                    % %Determine the midpoint value for each char
                    % Xs_bs_below_full = [Xs_bs_below, ones(params.Nz+1, 1)];
                    % Xs_bs_below_mid = (Xs_bs_below_full(:,1:end-1) + Xs_bs_below_full(:,2:end))./2; %Actual value of the cell 
                    % Xs_bs_below_layer_mid = Xs_bs_below_mid(zind+1,:);
                    % 
                    % 
                    % %Determine which bin to assign it ot
                    % inds_X_equal = Xs_bs_below_layer == X_in;
                    % [~,inds_X_closest] = min(abs(Xs_bs_below_layer - X_in));
                    % 
                    % 
                    % %Evaluate if there are multiple valid bins
                    % if length(find(inds_X_equal)) > 1 %Case where it is equal to multiple bins - use lowest conversion
                    %     minind = Xs_bs_order(min(find(inds_X_equal)));
                    %     inds_below(ii) = Xinds_below(minind);
                    % elseif ~all(inds_X_equal) & length(inds_X_closest) == 1 %Not equal to any - 
                    %     ind = Xs_bs_order(inds_X_closest);
                    %     inds_below(ii) = Xinds_below(ind);
                    % else
                    %     error('Uh Oh.');
                    % 
                    % end


                    %inds_below(ii) = find(params.zinds == zind_below & params.xinds == xind...
                    %& params.Tinds == Tinds_below(ii));
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
                Tcs_above = 1;

                %Convert into overall index (account for space/mass)
                zind_above = zind + 1;
                inds_above = zeros(size(Tinds_above));
                for ii = 1:length(Tinds_above)


                    %Identify the bins above this one
                    XTinds_above = find(params.zinds == zind_above & params.xinds == xind);

                    %Identify the conversion bins above this one
                    Xinds_above = find(params.zinds == zind_above & params.xinds == xind...
                        & params.Tinds == Tinds_above(ii));

                    %Pull boundary values
                    Tc_above = zeros(size(Xinds_above));
                    Tb_above = Tc_above;
                    Xc_above = Tc_above;
                    Xb_above = Tc_above;
                    Xb_final = Tc_above;
                    for ia = 1:length(Xinds_above)

                        %Pull indices
                        XTind = Xinds_above(ia);
                        ix = params.Xinds(XTind);
                        it = params.Tinds(XTind);

                        %Determine the temperatures above
                        Tcs_above = params.T_z.chars.merged{xind}{ix}.Tcz(:,zind_above);
                        active_inds = Tcs_above > 0;
                        Tcs_above = Tcs_above(active_inds);
                        Tc_above(ia) = Tcs_above(it);

                        %Determine the conversions above
                        Xcs_above = params.T_z.chars.merged{xind}{ix}.Xcz(:,zind_above);
                        Xcs_above = Xcs_above(active_inds);
                        Xc_above(ia) = Xcs_above(it);
                        Xbs_above = params.T_z.chars.merged{xind}{ix}.Xbz(:,zind_above);
                        Xbs_above = Xbs_above(active_inds);
                        Xb_above(ia) = Xbs_above(it);
                        Xbs_final = params.T_z.chars.merged{xind}{ix}.Xbz(:,end);
                        Xbs_final = Xbs_final(active_inds);
                        Xb_final(ia) = Xbs_final(it);

                    end

                    %Determine bin above 
                    if length(Xinds_above) == params.chars.N_Xis %Simple case - direct mapping
                        inds_above(ii) = find(params.zinds == zind_above & params.xinds == xind...
                        & params.Tinds == Tinds_above(ii) & params.Xinds == Xind);
                    else
                        Xinds_X_above = params.Xinds(Xinds_above);
                        inds_X_equal = Xb_above == X_out;
                        [~,inds_X_closest] = min(abs(Xb_above - X_out));
                        if any(Xinds_X_above == Xind & inds_X_equal) 
                            inds_same = Xinds_above(Xinds_X_above == Xind);
                            inds_above(ii) = inds_same;
                        elseif all(Xb_above < X_out)
                            inds_above(ii) = Xinds_above(end);
                        elseif all(Xb_above > X_out)
                            inds_above(ii) = Xinds_above(1);
                        else
                            inds_above(ii) = Xinds_above(inds_X_closest);
                        end
                       
                    end

                        %

    
    %                     %Pull Conversion values at boundaries
    %                     Xs_bs_above = zeros(params.Nz+1, params.chars.N_Xis);
    %                     for ix = 1:params.chars.N_Xis
    %                         Xs_bs_above(:,ix) = params.T_z.chars.Xs_bs{xind}{Tinds_above(ii)}{ix};
    %                     end
    %                     [Xs_bs_f, Xs_bs_order] = sort(Xs_bs_above(end,:));
    %                     Xs_bs_above_orig = Xs_bs_above;
    %                     Xs_bs_above = Xs_bs_above(:, Xs_bs_order);
    %                     Xs_bs_above_layer = Xs_bs_above(zind+1,:);
    % 
    % 
    %                     %Determine the midpoint value for each char
    %                     Xs_bs_above_full = [Xs_bs_above, ones(params.Nz+1, 1)];
    %                     Xs_bs_above_mid = (Xs_bs_above_full(:,1:end-1) + Xs_bs_above_full(:,2:end))./2; %Actual value of the cell 
    %                     Xs_bs_above_layer_mid = Xs_bs_above_mid(zind+1,:);
    % 
    %                     %Determine which above bin to map this bin to
    %                     %min_X_above = min(Xs_bs_above_layer(Xs_bs_above_layer >= X_out));
    %                     %max_X_below = max(Xs_bs_above_layer(Xs_bs_above_layer <= X_out));
    %                     %inds_X_above = Xs_bs_above_layer == min_X_above;
    %                     %inds_X_below = Xs_bs_above_layer == max_X_below;
    %                     inds_X_equal = Xs_bs_above_layer == X_out;
    %                     [~,inds_X_closest] = min(abs(Xs_bs_above_layer - X_out));
    % 
    %                     %inds_X_above = Xs_bs_above_layer >= X_out;
    %                     %inds_X_below = Xs_bs_above_layer <= X_out;
    % 
    %                     %Evaluate if there are multiple valid bins
    %                     %if all(inds_X_equal)
    % %
    %                     %end
    % 
    %                     if length(find(inds_X_equal)) > 1 %Case where it is equal to multiple bins - use lowest conversion
    %                         minind = Xs_bs_order(min(find(inds_X_equal)));
    %                         inds_above(ii) = Xinds_above(minind);
    %                     elseif ~all(inds_X_equal) & length(inds_X_closest) == 1 %Not equal to any - 
    %                         ind = Xs_bs_order(inds_X_closest);
    %                         inds_above(ii) = Xinds_above(ind);
    %                     else
    %                         error('Uh Oh.');
    % 
    %                     end
    % 
    %                 end
    % 
    % 
    % 
    % 
    %                 % inds_X_valid = find(inds_X_above & inds_X_below);
    %                 % 
    %                 % %Map to the bin that gives the lowest conversion to
    %                 % %be most conservative
                    % if isempty(inds_X_valid) %No bin it fits in
                    %     if all(inds_X_below) %Assign to bin that results in highest conversion
                    %         [~, maxind] = max(Xs_bs_above_layer(end,:));
                    %         inds_above(ii) = Xinds_above(maxind);
                    %     elseif all(inds_X_above) %Assign to bin that results in lowest conversion
                    %         [~, minind] = min(Xs_bs_above_layer(end,:));
                    %         inds_above(ii) = Xinds_above(minind);
                    %     else
                    %         error('No valid bin. Programming issue.')
                    %     end
                    % else
                    %     Xs_f_valid = Xs_bs_above(end,inds_X_valid);
                    %     minind = find(Xs_f_valid == min(Xs_f_valid));
                    %     inds_above(ii) = Xinds_above(minind);
                    % end



          
                end
                
                %Store result
                params.inds_above{i} = inds_above;
            end



        elseif params.heat.active

            %Calculate counts
            NT = unique(params.Tinds(params.xinds == xind & params.zinds == zind));
            NX = 1; %TO BE UPDATED
    
            %Pull temperature characteristics
            merged = params.T_z.chars.merged{xind}{1};

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