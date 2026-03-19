function T_bar = CalcTbar(y_t, y_t_m, iz, im, params)

%Run without input
if nargin < 1
    PBM_postprocess();
end



%Calculations


    if params.heat.active %Weighted average

        %Pull overall cell value
        ind_m = find(params.zinds_m(:) == iz & params.xinds_m(:) == im);
        N = y_t_m(ind_m);

        %Pull merged values
        merged = params.T_z.chars.merged{im};

        %Pull corresponding indices
        subinds = find(params.xinds == im & params.zinds == iz);
        subys = y_t(subinds);

        %Identify curve of slowest mass
        Ts_slowest = params.T_z.chars.Ts_cs{end}{1}{1};

        %Pull trailing inds and identify inactive ones
        NTs = length(params.T_z.chars.merged{im}.bins_in(:,1));
        %inds_trail = params.T_z.chars.itrails{im};
        Tczs = params.T_z.chars.merged{im}.Tcz(:,iz); Tczs = flipud(Tczs);
        inactive_inds = Tczs < (Ts_slowest(iz) + 3) & Tczs > 0;
        inactive_inds = find(inactive_inds);
        %inactive_inds = NTs - inactive_inds;
        
        %Identify active bins
        if iz == 1
            active_inds = params.T_z.chars.merged{im}.bins_in(:,1);
            %active_inds = flipud(active_inds);
        else
            bins_in = params.T_z.chars.merged{im}.bins_in(:,iz);
            bins_in = flipud(bins_in);
            active_inds = Tczs > -1; active_inds = active_inds(1:end-1);
            trail_inds = Tczs < (Ts_slowest(iz) + 3) & Tczs > 0;
            trail_inds = trail_inds(2:end);
            active_inds = active_inds & ~trail_inds; 
            active_inds = bins_in(active_inds);

            active_inds = unique(bins_in);

        end
        %active_inds = active_inds(active_inds > max(inds_trail));
        %active_inds = active_inds - 1;

        %Calculate fraction for each temperature
        if N > 0

            yfrac = subys(active_inds)./N;
        else
            yfrac = zeros(size(active_inds));
        end


        %Pull corresponding temperatures
        Tcz = params.T_z.chars.merged{im}.Tcz;
        if iz == 1
            Ts = zeros(1, length(active_inds)); %length(params.T_z.chars.Ts_cs{im}));
            ix = 1; %UPDATE UPDATE UPDATE 
            for i = 1:length(Ts)
                ind = active_inds(i);
                Ts(i) = params.T_z.chars.Ts_cs{im}{ind}{ix}(iz);
            end
            Ts = fliplr(Ts);
        else
            Ts = Tcz(find(Tcz(:,iz) > 0), iz);
            Ts = (Ts(1:end-1) + Ts(2:end))./2;        
            Ts = flipud(Ts);%Flip Ts
        end

        

        %Calculate weighted mean
        if length(Ts) == length(yfrac)
            T_bar = sum(Ts(:) .* yfrac(:));
        else
            error('Ts and yfrac length do not match. Programming error.');
        end
        


        
        x = 1;


    else %Assume uniform distribution of temperatures and conversions

        try

            %Pull appropriate temperature estimates
            Ts_cs = params.T_z.chars.Ts_cs{im};

            %Create a storage matrix
            T_mat = zeros(length(Ts_cs), params.chars.N_Xis);

            for it = 1:length(Ts_cs)
                Ts_c = Ts_cs{it};
                for ix = 1:params.chars.N_Xis
                    T_mat(it, ix) = Ts_c{ix}(iz);
        
                end
            end
    
            T_bar = mean(mean(T_mat));
        catch
            T_bar = params.T_liq;
        end

    end
    

end