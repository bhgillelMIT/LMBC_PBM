%This function takes in a set of source terms calculated for the mean
% conditions of each representative mass (im = 1:params.Nms) in each
% spatial cell (iz = 1:params.Nz). It then distributes the source terms
% between the subcells (temperature and conversion) based on their
% distribution in the source representative mass. The function must 
% an energy and mass balance to keep the system physical.

function h = DistSourceTerms(y, h_m, cadd, csub, badd, bsub, params)

    %Settings
    reltol = 1E-6;
    debug = true;
    if debug
        %debugfig = figure();
    end

    %Debug message and counter
    if params.debug
        t_start = cputime;
        fprintf('--Distributing Source Terms Start\n');
    end

    %Initialize output
    h = params.h;

    %Extract relevant variables
    xinds_m = params.xinds_m;
    zinds_m = params.zinds_m;
    xinds = params.xinds;
    zinds = params.zinds;
    Tinds = params.Tinds;
    Xinds = params.Xinds;
    cmats = params.cmats;
    bmats = params.bmats;
    
    switch params.heat.source_mode
        case 'adaptive'
   
    
        %Iterate through each spatial cell 
        i = 1; % # - iteration counter
        for iz = 1:params.Nz

            %Debug statement
            if params.debug
                fprintf('---Layer %d/%d Start:\n', iz, params.Nz)
            end

            %Pull relevant Ts
            zinds_z = zinds == iz;
            relinds_z = find(zinds_z); %indices in this spatial cell
            relTs_z = params.Ts(zinds_z);
            relits_z = params.Tinds(zinds_z);
            relXs_z = params.Xs(zinds_z);
            relixs_z = params.Xinds(zinds_z);
            relys = y(relinds_z);

            %Calcualte largest T and X
            zm = params.zms(iz);
            T_largest = interp1(params.T_z.chars.largest.zs, params.T_z.chars.largest.Ts, zm);
            X_largest = interp1(params.T_z.chars.largest.zs, params.T_z.chars.largest.Xs, zm);


            %Iterate through each mass - receiving influx (sink)
            for im = 1:params.Nms

                

                %Define overall index
                xind = (iz-1)*params.Nms + im; 

                %Identify indices
                subind_m = find(xinds_m == im & zinds_m == iz);
                subinds = find(xinds == im & zinds == iz);
                

                %Pull source terms specific to this mass
                h_im = h_m(xind);
                cadd_m = cadd(xind);
                csub_m = csub(xind);
                badd_m = badd(xind);
                bsub_m = bsub(xind);
                hadd_m = cadd_m + badd_m;
                hsub_m = csub_m + bsub_m;

                %Pull values
                snkN = params.Ns_m(subind_m);
                snkNs = y(subinds, :); %All subbins for this mass
                snkTs = params.Ts(subinds);
                snkXs = params.Xs(subinds);
                snkit = params.Tinds(subinds);
                snkix = params.Xinds(subinds);
                snkTinds = params.Tinds(subinds);
                snkXinds = params.Xinds(subinds);
                snkT_aboves = params.Ts_above(subinds);
                snkX_aboves = params.Xs_above(subinds);
                snkT_finals = params.Ts_final(subinds);
                snkX_finals = params.Xs_final(subinds);

                %Determine which of the bins are trialing 
                snk_trailing_T = abs(snkTs - T_largest)./T_largest < reltol;
                snk_trailing_X = abs(snkXs - X_largest)./X_largest < reltol;
                snk_trailing = snk_trailing_T & snk_trailing_X;

                %Pull the active matrix for the cells
                snk_actives = cell2mat(params.cell_active(subinds));

                %Identify the trailing bin to assign other trailing points
                %to
                [X_trailing_max, maxind] = max(snkX_finals(snk_trailing));
                snk_trailing_inds = find(snk_trailing);
                snk_trailing_ind = snk_trailing_inds(maxind); %Bin to assign incoming trailing points to

                %Identify smallest and largest bins - the one with the
                %lowest final conversion
                if any(snk_trailing)
                    [X_trailing_min, minind] = min(snkX_finals(snk_trailing));
                    snk_min_inds = find(snkX_finals == X_trailing_min);
                    snk_min_ind = snk_min_inds(1);
                else
                    snk_min_ind = 1;
                end

                [X_max, maxind] = max(snkX_finals(~snk_trailing));
                snk_max_inds = find(snkX_finals == X_max);
                snk_max_ind = snk_max_inds(end);


                %Calculate fluxes
                h_iT_in = zeros(size(snkTs));
                h_iT_out = h_iT_in;
                %h_iT_in = zeros(size(snkNs)); %Vector keeping track of flux into each cell 

                %hadd = params.badd(xind) + params.cadd(xind); %fluxes in for this mass
                %hsub = params.bsub(xind) + params.csub(xind); %fluxes out for this mass

               
                %Distribute incoming source - only needed if cadd_m or
                %badd_m is greater than 0. Otherwise it is just a flux out
                if hadd_m > 0

                    %Load breakage and coalescence fluxes
                    if params.coalesce.active && params.break.active
                        cmat = cmats{iz}; cvec = cmat.src(:,im);
                        bmat = bmats{iz}; bvec = bmat(:,im); 
                    elseif params.coalesce.active
                        cmat = cmats{iz}; cvec = cmat.src(:,im);
                        bmat = zeros(params.Nms); bvec = bmat(:,im);
                    elseif params.break.active
                        bmat = bmats{iz}; bvec = bmat(:,im);
                        cmat = zeros(params.Nms); cvec = cmat.src(:,im);
                    else
                        bmat = zeros(params.Nms); bvec = bmat(:, im);
                        cmat = zeros(params.Nms); cvec = cmat.src(:, im);
                    end
    
                    %Check conservation
                    c_src = sum(cvec);
                    c_snk = sum(cmat.src(im,:));
                    b_src = sum(bvec);
                    b_snk = sum(bmat(im,:));
    
                   
    
                    %Calculate normal vectors
                    if any(bvec > 0)
                        bvec_norm = bvec./sum(bvec);
                    else
                        bvec_norm = ones(size(bvec));
                    end
                    if any(cvec > 0)
                        cvec_norm = cvec./sum(cvec);
                    else
                        cvec_norm = zeros(size(cvec));
                    end

                    %Convert to mass flux 



                    %Identify masses contributing to this mass
                    srcinds_m = find(cvec > 0 | bvec > 0); %Representative mass indices which have fluxes into this representative mass (im)
                    
                    %Iterate through other sizes - providing influx
                    %(source)
                    for i = 1:length(srcinds_m)

                        %Pull index and normalized value
                        ind = srcinds_m(i); %Mass bin of this source
                        bnorm = bvec_norm(ind);
                        cnorm = cvec_norm(ind);
                        babs = bvec(ind); %Absolute flux from this 
                        cabs = cvec(ind);
                        sabs = babs + cabs; %Total influx of source terms
        
                        %Pull temperatures and conversions
                        srcinds_T = find(xinds(relinds_z) == ind); %Indices of the numeric density corresponding to the source mass
                        srcTs = relTs_z(srcinds_T);
                        srcit = relits_z(srcinds_T);
                        srcTs_above = params.Ts_above(relinds_z(srcinds_T));
                        srcXs = relXs_z(srcinds_T);
                        srcix = relixs_z(srcinds_T);
                        srcXs_above = params.Xs_above(relinds_z(srcinds_T));

                        %Determine the actual bins the 
                        srcNs = relys(srcinds_T);
                        if all(srcNs == 0)
                            srcNfracs = srcNs./(sum(srcNs)+1E-32);
                        else
                            srcNfracs = srcNs./(sum(srcNs));
                        end

                        % %Debug plot
                        % if debug
                        %     figure(debugfig)                            
                        %     plot(srcTs, srcNs, 'o-', 'LineWidth', 2); hold on;
                        %     xlabel('Temperature (K)'); ylabel('Numeric Density'); 
                        %     title('Source Distribution'); 
                        %     grid on; grid minor; axis square;
                        %     cla;
                        % 
                        % end

                        %Iterate through each source 
                        for it_src = 1:length(srcTs)
                            T = srcTs(it_src);
                            X = srcXs(it_src);
                            N = srcNs(it_src);
                            Nfrac = srcNfracs(it_src);

                            %Identify the adjacent cells in the sink mass
                            if Nfrac > 0
                                if params.react.active

                                    %Categorize the point
                                    trailing = abs(T-T_largest) < 0.0001 & abs(X-X_largest) < 0.0001; %logical - is the source bin a trailing bin
                                    trailing = trailing | (T <= T_largest & X <= X_largest);
                                    smaller_T = all(T < snkTs);
                                    smaller_X = all(X < snkXs);
                                    smaller = smaller_T & smaller_X;
                                    larger_T = all(T > snkTs);
                                    larger_X = all(X > snkXs);
                                    larger = larger_T & larger_X;
                                    matchingT = abs(T - snkTs)/T < reltol;
                                    matchingX = abs(X - snkXs)/X < reltol;
                                    matchingBOTH = matchingT & matchingX;

                                    %Distribute 
                                    if trailing %Distribute to the bin that reaches the largest final conversion
                                        h_iT_in(snk_trailing_ind) = h_iT_in(snk_trailing_ind) + sabs * Nfrac;
                                    else %Determine which bin to assign it to
                                        if smaller
                                            h_iT_in(snk_min_ind) = h_iT_in(snk_min_ind) + sabs * Nfrac;
                                        elseif larger
                                            h_iT_in(snk_max_ind) = h_iT_in(snk_max_ind) + sabs * Nfrac;
                                        else %Determine which bin it fits into of the active bins
                        
                                            

                                            %Consider only the active bins
                                            actives = snk_actives(:,1);
                                            actives_inds = find(actives);
                                            
                                            snkXs_active = snkXs(actives);
                                            snkTs_active = snkTs(actives);

                                            %Determine if there is one bin
                                            %that matches perfectly
                                            X_match = snkXs_active == X;
                                            T_match = snkTs_active == T;
                                            BOTH_match = X_match & T_match;
                                            
                                            N_match = length(find(BOTH_match));
                                            if N_match == 1
                                                ind = actives_inds(BOTH_match);
                                                h_iT_in(ind) = h_iT_in(ind) + sabs * Nfrac;
                                            elseif N_match > 1

                                                %Identify the indices that
                                                %result in the highest
                                                %conversion
                                                X_finals = snkX_finals(actives);
                                                X_finals_BOTH = X_finals(BOTH_match);
                                                X_finals_max = max(X_finals_BOTH);
                                                active_inds = find(actives);
                                                BOTH_match_inds = find(BOTH_match);
                                                max_inds = X_finals_BOTH == X_finals_max;
                                                max_inds = find(X_finals_BOTH(find(max_inds)));
                                                snk_inds = active_inds(BOTH_match_inds(max_inds));
                                                N_inds = length(snk_inds);


                                                %Distribute between them
                                                %equally
                                                dist_frac = 1/N_inds;
                                                for is = 1:N_inds
                                                    snk_ind = snk_inds(is);
                                                    h_iT_in(snk_ind) = h_iT_in(snk_ind) + dist_frac*sabs*Nfrac;
                                                end


                                                x = 1;
                                            else %Identify the bin it is closest to

                                                err_X = cell(1, params.chars.N_Xis);
                                                err_T = cell(1, params.chars.N_Xis);
                                                err_T_diff = cell(1, params.chars.N_Xis);
                                                err_min = zeros(1, params.chars.N_Xis);
                                                err_X_diff = cell(1, params.chars.N_Xis);
                                                err_total = cell(1, params.chars.N_Xis);
                                                for ix = 1:params.chars.N_Xis

                                                    %Define sub indices 
                                                    subinds_X = snkXinds == ix;

                                                    %Pull Ts
                                                    Ts_X = snkTs(subinds_X);
                                                    Xs_X = snkXs(subinds_X);

                                                    %Calculate err
                                                    err_X{ix} = X - Xs_X;
                                                    err_X_diff{ix} = diff(err_X{ix});
                                                    err_T{ix} = T - Ts_X;
                                                    err_T_diff{ix} = diff(err_T{ix});
                                                    err_total{ix} = sqrt((err_X{ix}./X).^2 + (err_T{ix}./T).^2);
                                                    err_min(ix) = min(err_total{ix});

                                                end
            
                                                %Identify the actual bins
                                                %to divide between
                                                [min_err, ix_min] = min(err_min);
                                                err_totals = err_total{ix_min};
                                                [min_err, iT_min] =  min(err_totals);
                                                if iT_min == 1
                                                    snk_lo = find(iT_min == snkTinds & ix_min == snkXinds);
                                                    snk_hi = find(iT_min+1 == snkTinds & ix_min == snkXinds);
                                                elseif iT_min == length(err_totals)
                                                    snk_lo = find(iT_min-1 == snkTinds & ix_min == snkXinds);
                                                    snk_hi = find(iT_min == snkTinds & ix_min == snkXinds);
                                                else
                                                    err_neg = err_totals(iT_min-1);
                                                    err_pos = err_totals(iT_min+1);
                                                    if err_neg < err_pos
                                                        snk_lo = find(iT_min-1 == snkTinds & ix_min == snkXinds);
                                                        snk_hi = find(iT_min == snkTinds & ix_min == snkXinds);
                                                    else
                                                        snk_lo = find(iT_min == snkTinds & ix_min == snkXinds);
                                                        snk_hi = find(iT_min+1 == snkTinds & ix_min == snkXinds);
                                                    end
                                                end
                                                
                                                %Distribute
                                                snk_ind = find(iT_min == snkTinds & ix_min == snkXinds);
                                                h_iT_in(snk_ind) = h_iT_in(snk_ind) + sabs*Nfrac;


                
                                                % 
                                                % %Identify T bounds
                                                % T_lower = max(snkTs_active(snkTs_active <= T));
                                                % T_upper = min(snkTs_active(snkTs_active >= T));
                                                % T_lower_inds = find(snkTs_active == T_lower);
                                                % T_upper_inds = find(snkTs_active == T_upper);
                                                % 
                                                % 
                                                % %Identify X bounds
                                                % X_lower = max(snkXs_active(snkXs_active <= X));
                                                % X_upper = min(snkXs_active(snkXs_active >= X));
                                                % X_lower_inds = find(snkXs_active == X_lower);
                                                % X_upper_inds = find(snkXs_active == X_upper);

                                                % %Identify closest two bins
                                                % snk_lo = max(find((snkTs(actives) - T) <= 0));
                                                % snk_hi = min(find((snkTs(actives) - T) > 0));
                                                % 
                                                % 
                                                % %Calculate division using 
                                                % eta_lo = (snkTs(snk_hi) - T)./(snkTs(snk_hi) - snkTs(snk_lo));
                                                % eta_hi = 1 - eta_lo;
                                                % 
                                                % %Distribute
                                                % h_iT_in(snk_lo) = h_iT_in(snk_lo) + eta_lo * sabs * Nfrac;
                                                % h_iT_in(snk_hi) = h_iT_in(snk_hi) + eta_hi * sabs * Nfrac;


                                            
                                            end

                                            


        
                                            
    
                                    

                                            x = 1;

                                        end

                                    end

                                    % 
                                    % %Determine which snk bins are active or
                                    % %activating
                                    % active_now = snk_actives(:,1);
                                    % if iz ~= params.Nz
                                    %     active_next = snk_actives(:,2);
                                    % else
                                    %     active_next = false(size(active_now));
                                    % end
                                    % activating = active_next & ~active_now;
                                    % 
                                    % matching_active = matchingBOTH & active_now;
                                    % matching_activating = matchingBOTH & activating;

                                    
                                    
    
                                    % if T < min(snkTs) %Case where it goes into the last trailing bin
                                    % 
                                    %     h_iT_in(1) = h_iT_in(1) + sabs * Nfrac;
                                    % 
                                    % elseif T >= max(snkTs)
                                    % 
                                    %     h_iT_in(end) = h_iT_in(end) + sabs * Nfrac;
                                    % 
                                    % elseif any(matchingBOTH) %Indicates likely a trailing bin
                                    % 
                                    %     %Determien which sink is largest at the
                                    %     %end
                                    %     if length(find(matchingBOTH)) > 1
                                    % 
                                    % 
                                    %         Xs_final_valid = cell2mat(snkX_aboves(matching_activating));
                                    %         [largest_X, largest_ind] = max(Xs_final_valid(:,end));
                                    % 
                                    %         %Determine if there remain multiple -
                                    %         %ONLY DO THIS IF TRAILING
                                    %         %bins or if there are too few
                                    %         if length(largest_ind) == 0
                                    %             error('No valid index to assign to. Programming error.');
                                    %         elseif length(largest_ind) > 1
                                    %             assign_ind = find(matching_activating);
                                    %             assign_ind = assign_ind(largest_ind(1)); %Just take the first one for simplicity
                                    %         else
                                    %             assign_ind = find(matching_activating);
                                    %             assign_ind = assign_ind(largest_ind);
                                    %         end
                                    % 
                                    % 
                                    %         %Identify the next X bin
                                    % 
                                    %     elseif length(find(matchingBOTH)) == 0
                                    %         error('Invalid option.')
                                    %     else
                                    %         assign_ind = find(matchingBOTH);
                                    %     end
                                    % 
                                    % 
                                    % 
                                    %     x = 1;
        
                                    % else %Identify closest two bins and distribute between them
                                    % 
                                    %     %Identify closest two bins
                                    %     snk_lo = max(find((snkTs - T) <= 0));
                                    %     snk_hi = min(find((snkTs - T) > 0));
                                    % 
                                    %     %Calculate division using 
                                    %     eta_lo = (snkTs(snk_hi) - T)./(snkTs(snk_hi) - snkTs(snk_lo));
                                    %     eta_hi = 1 - eta_lo;
                                    % 
                                    %     %Distribute
                                    %     h_iT_in(snk_lo) = h_iT_in(snk_lo) + eta_lo * sabs * Nfrac;
                                    %     h_iT_in(snk_hi) = h_iT_in(snk_hi) + eta_hi * sabs * Nfrac;
                                    % 
                                    % end
    
    
                                else
        
                                    if T < min(snkTs) %Case where it goes into the last trailing bin
                                      
                                        h_iT_in(1) = h_iT_in(1) + sabs * Nfrac;
                                        
                                        %Adjust temperature to preserve energy -
                                        %ADD
                                   
        
                                    elseif T >= max(snkTs)
        
                                        h_iT_in(end) = h_iT_in(end) + sabs * Nfrac;
        
                                    else %Identify closest two bins and distribute between them
        
                                        %Identify closest two bins
                                        snk_lo = max(find((snkTs - T) <= 0));
                                        snk_hi = min(find((snkTs - T) > 0));
        
                                        %Calculate division using 
                                        eta_lo = (snkTs(snk_hi) - T)./(snkTs(snk_hi) - snkTs(snk_lo));
                                        eta_hi = 1 - eta_lo;
        
                                        %Distribute
                                        h_iT_in(snk_lo) = h_iT_in(snk_lo) + eta_lo * sabs * Nfrac;
                                        h_iT_in(snk_hi) = h_iT_in(snk_hi) + eta_hi * sabs * Nfrac;
        
                                        x = 1;
                                        
        
                                    end
                                end
                            end

                            %Detect NaNs
                            if any(isnan(h_iT_in))
                                warning('NaNs');
                            end



                        end

                        %Store results
                        


                        %Calculate flux in
                        x = 1;
        
                    end 

                else
                    x = 1;
                end

                %Calculate the flux out
                if all(snkNs == 0)
                    h_iT_out = (csub_m + bsub_m) * (snkNs./(sum(snkNs)+1E-16));
                else
                    h_iT_out = (csub_m + bsub_m) * (snkNs./(sum(snkNs)));
                end


                %Concatenate flux vector to 
                h(subinds) = h_iT_in - h_iT_out;

                %Check error
                if debug

                    %Calculate absolute and relative error
                    hadd_T = sum(h_iT_in);
                    flux_aerr = abs(sum(h_iT_in) - hadd_m);
                    flux_rerr = flux_aerr/(hadd_m+1E-16);
                    if flux_rerr > 0.0001
                        warning('Issue with DistSourceTerms. Flux doesn''t balance');
                    end

                    %Debug statement
                    fprintf('----Mass (im = %d) Distribution (h_m = %E; h_T = %E; rerr = %E)\n', im, hadd_m, hadd_T, flux_rerr);

                end



                
                % if hi > 1
                %     x = 1;
                % end

                x = 1;


                %Update iteration counters
                i = i + 1;

            end

            %Debug check
            if debug

                % %Check sum of source terms
                % h_total = sum(h);
                % h_m_total = sum(h_m);
                % h_aerr = abs(h_total - h_m_total);
                % h_rerr = h_aerr/(h_m_total+1E-16);
                % fprintf('---Layer %d/%d - Source Term Check: Total h = %0.4e, Total h_m = %0.4e, Abs Error = %0.4e, Rel Error = %0.4e\n', iz, params.Nz, h_total, h_m_total, h_aerr, h_rerr);

                %Plot the distributions for comparison
                % figure();
                % subplot(1,2,1);
                % plot(h, 'k-'); hold on;
                % 
                % subplot(1,2,2);
                % plot(h_m, 'b-'); hold on;

            end



        end

        case 'simple'

            %Iterate through each size/spatial cell, pull the net flux into that cell, and distribute it based on the original breakdown for the temperature distributions
            for iz = 1:params.Nz
                    for im = 1:params.Nms
    
                        %Define overall index
                        xind = (iz-1)*params.Nz + im; 
    
                        %Identify indices
                        subind_m = find(xinds_m == im & zinds_m == iz);
                        subinds = find(xinds == im & zinds == iz);
    
                        %Pull values
                        subN = params.Ns_m(subind_m);
                        subNs = y(subinds); %All subbins for this mass
                        hi = h_m(xind);
                        hadd = params.badd(xind) + params.cadd(xind); %fluxes in for this mass
                        hsub = params.bsub(xind) + params.csub(xind); %fluxes out for this mass
    
                        if hi > 0
    
                            %Calculate distribution of source terms based on original distribution of temperatures for this mass
                            Tdist = subNs./sum(subNs(:));
                            h(subinds) = h(subinds) + hi.*Tdist;
    
                        end
    
                end
           end

        otherwise %Output error
            error('Invalid heat source mode specified. Options are simple or adaptive.');
    end


    %Debug message
    if params.debug
        t_end = cputime;
        t_req = t_end - t_start;
        fprintf('--Distributing Source Terms End (t = %0.4f s)\n', t_req);
    end




end