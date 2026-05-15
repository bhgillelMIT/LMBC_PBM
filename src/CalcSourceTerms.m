function [h, h_m, cadd, csub, badd, bsub, params] = CalcSourceTerms(t, F, params)

    %global params
    global output


    if t < params.sol.src_delay
        h = zeros(size(F));
        h_m = zeros(1, params.Nms * params.Nz);
        cadd = h_m; csub = h_m; badd = h_m; bsub = h_m;
    else
        params.src.its = params.src.its + 1;
        cadd = zeros(length(F), 1); csub = cadd;
        badd = cadd; bsub = cadd;
        if params.coalesce.active
            [cadd, csub, cmats] = Coalescence(params.Ns_m, params); %Coalescence(Fs, params );  
            cadd = cadd * params.coalesce.damper; %Damper is only used for debugging
            csub = csub * params.coalesce.damper;
            params.cmats = cmats; %Store for heat/reaction calculations
        end
        if params.break.active
            if params.sol.solve_details
                [badd, bsub, bmats] = Breakage(params.Ns_m, params); %Breakage(Fs, params);
                params.break.badd = badd;
                params.break.bsub = bsub;
                params.bmats = bmats;
            else
                badd = params.break.badd;
                bsub = params.break.bsub;
            end
            badd = params.break.damper * badd;
            bsub = params.break.damper * bsub;
        end

        %Distribute
        if params.heat.active
            h = params.h;

            h_m = cadd(:) - csub(:) + badd(:) - bsub(:);

            %Store in params
            params.cadd = cadd; params.csub = csub;
            params.badd = badd; params.bsub = bsub;

            %Function to distribute proportionally
            if any(abs(h_m) > 1E-12)
                h = DistSourceTerms(F, h_m, cadd, csub, badd, bsub, params);
            else
                h = params.h;
            end


            %Iterate through bins and allocate
            ind = 1; %index for coalescence and breakage vectors 
            for iz = 1:params.Nz
                for im = 1:params.Nms
                    subinds = find(params.zinds == iz & params.xinds == im);
                    h(subinds) = params.Ns_fracs(subinds) .* (cadd(ind) - csub(ind)...
                        + badd(ind) - bsub(ind));
                    ind = ind + 1;
                end
            end
        else
            h = cadd(:) - csub(:) + badd(:) - bsub(:);
            h_m = h;
        end

         
    end


    if any(isnan(h)) || any(~isreal(h))
        x = 1;
    end


end