%This function takes in a set of source terms calculated for the mean
% conditions of each representative mass (im = 1:params.Nms) in each
% spatial cell (iz = 1:params.Nz). It then distributes the source terms
% between the subcells (temperature and conversion) based on their
% distribution in the source representative mass. The function must 
% an energy and mass balance to keep the system physical.

function h = DistSourceTerms(y, h_m, cadd, csub, badd, bsub, params)

    %Settings
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
            relinds_z = find(zinds == iz); %indices in this spatial cell
            relTs_z = params.Ts(zinds == iz);
            relys = y(relinds_z);

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
                        cvec_norm = ones(size(cvec));
                    end

                    %Convert to mass flux 



                    %Identify masses contributing to this mass
                    srcinds_m = find(cvec > 0 | bvec > 0); %Representative mass indices which have fluxes into this representative mass (im)
                    
                    %Iterate through other sizes - providing influx
                    %(source)
                    for i = 1:length(srcinds_m)

                        %Pull index and normalized value
                        ind = srcinds_m(i);
                        bnorm = bvec_norm(ind);
                        cnorm = cvec_norm(ind);
                        babs = bvec(ind); %Absolute flux from this 
                        cabs = cvec(ind);
                        sabs = babs + cabs; %Total influx of source terms
        
                        %Pull temperatures
                        srcinds_T = find(xinds(relinds_z) == ind); %Indices of the numeric density corresponding to the source mass
                        srcTs = relTs_z(srcinds_T);
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
                            N = srcNs(it_src);
                            Nfrac = srcNfracs(it_src);

                            %Identify the adjacent cells in the sink mass
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