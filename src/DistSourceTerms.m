%This function takes in a set of source terms calculated for the mean
% conditions of each representative mass (im = 1:params.Nms) in each
% spatial cell (iz = 1:params.Nz). It then distributes the source terms
% between the subcells (temperature and conversion) based on their
% distribution in the source representative mass. 

function h = DistSourceTerms(y, h_m, params)

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
    
    if strcmp(params.heat.source_mode, 'adaptive')
    
        %Iterate through each size
        i = 1;
        for iz = 1:params.Nz

            %Pull relevant Ts
            relinds_z = find(zinds == iz);
            relTs_z = params.Ts(zinds == iz);

            for im = 1:params.Nms

                %Define overall index
                xind = (iz-1)*params.Nz + im; 

                %Identify indices
                subind_m = find(xinds_m == im & zinds_m == iz);
                subinds = find(xinds == im & zinds == iz);

                %Pull values
                subN = params.Ns_m(subind_m);
                subNs = y(subinds, :); %All subbins for this mass
                hi = h_m(xind);
                hadd = params.badd(xind) + params.cadd(xind); %fluxes in for this mass
                hsub = params.bsub(xind) + params.csub(xind); %fluxes out for this mass

                %Load breakage and coalescence fluxes
                if params.coalesce.active && params.break.active
                    cmat = cmats{iz}; cvec = cmat(:,im);
                    bmat = bmats{iz}; bvec = bmat(:,im); 
                elseif params.coalesce.active
                    cmat = cmats{iz}; cvec = cmat(:,im);
                    bmat = zeros(params.Nms); bvec = bmat(:,im);
                elseif params.break.active
                    bmat = bmats{iz}; bvec = bmat(:,im);
                    cmat = zeros(params.Nms); cvec = cmat(:,im);
                else
                    bmat = zeros(params.Nms); bvec = bmat(:, im);
                    cmat = zeros(params.Nms); cvec = cmat(:, im);
                end

                %Check conservation
                c_src = sum(cvec);
                c_snk = sum(cmat(im,:));
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

                if hi > 0

                    %Identify masses contributing to this mass
                    relinds = find(cvec > 0 | bvec > 0);
                    
                    %Iterate through other sizes
                    for i = 1:length(relinds)

                        %Pull index and normalized value
                        ind = relinds(i);
                        bnorm = bvec_norm(ind);
                        cnorm = cvec_norm(ind);
        
                        %Pull temperatures
                        relinds = find(xinds(relinds_z) == ind);
                        relTs = relTs_z(relinds);
                        relys = subys(relinds);

                        %Calculate flux in
                        x = 1
        
                    end 

                else
                    x = 1;
                end


                
                % if hi > 1
                %     x = 1;
                % end

                x = 1;


                %Update iteration counters
                i = i + 1;

            end
        end

    elseif strcmp(params.heat.source_mode, 'simple')

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
                subNs = y(subinds, :); %All subbins for this mass
                hi = h_m(xind);
                hadd = params.badd(xind) + params.cadd(xind); %fluxes in for this mass
                hsub = params.bsub(xind) + params.csub(xind); %fluxes out for this mass

                if hi > 0

                    %Calculate distribution of source terms based on original distribution of temperatures for this mass
                    Tdist = sum(subNs(:, Tinds), 2)./sum(subNs(:));
                    h(subinds) = h(subinds) + hi.*Tdist;

                end

            end
       end

     else
        error('Invalid heat source mode specified. Options are simple or adaptive.');
    end



end