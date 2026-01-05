function merged = MergeBins(params)

    %Plot settings
    ms = 18;
    fs = 18;
    lw = 2;

    %Debug settings
    Nm_max = params.disc.Nms;


    %Pull thresholds
    dT_cutoff = params.chars.dT_cutoff;


     %Create output cell array
    merged = cell(1,params.Nms);

    %Define zs
    zs_all = linspace(0, params.reactor.H, length(params.T_z.chars.Ts_zs{18}{1}{1}));

    %Iterate through mass
    

    for im = 1:Nm_max 

        %Pull characteristics
        Tzs = params.T_z.chars.Ts_zs{im};
        Tcs = params.T_z.chars.Ts_cs{im};

        %Determine number
        Nx = length(Tcs{1});

        %Define boundaries
        %bins_red = params.T_liq:-dT_cutoff: 


        % %Debug plot
        % if params.debug
        %     for ix = 1:Nx
        % 
        %     end
        % 
        % 
        % end


        
        %Iterate through initial conversion cases
        for ix = 1:1 %1:Nx

            %Debug plot
            if params.debug
                debug_fig = figure();
                subplot(1,2,1);
                for it = 1:length(Tcs)
                    plot(Tzs{it}{ix}, linspace(0, params.reactor.H, length(Tzs{it}{ix})), 'r-'); hold on;
                    plot(Tcs{it}{ix}, params.mesh.volcell_cents(:,2), 'r.', 'MarkerSize', ms); hold on;
                end
                title_str = sprintf('Nm = %d; Nx = %d;', im, ix);
                title(title_str);
            end


            %Pull current bin temperatures
            Nt = length(Tcs);
            Tcz = zeros(Nt, params.Nz);
            for it = 1:Nt
                Tcz(it,:) = Tcs{it}{ix}';
            end
            

            %Append liquid temperature to Tcz
            for iz = 1:params.Nz
                Tcz(Nt+1, iz) = params.T_Lz(params.mesh.volcell_cents(iz, 2));
            end
            Nt = Nt + 1;
            Tcz_total = Tcz; %Copy of the original bins

            %Allocate storage for new bounds
            inds_bin = 1:Nt;
            N_bins = zeros(1, params.Nz);

            %Iterate through Z-levels 
            inds_remaining = inds_bin;
            for iz = 1:params.Nz

                %Isolate this levels temperature bounds and count bins
                zc = params.zms(iz);
                Tc = Tcz(:,iz);
                %Nt = Nt + 1;
                N_bins_i = length(Tcz(:,iz))-1;

                %Determine initial bounds
                T_trail = Tcz(1,iz);
                T_liq_local = params.Tsz(iz);
                Tc_new = T_liq_local.*ones(1, params.Nz);

                %Determine how bins are trailing points
                inds_trail = find(Tcz_total(:,iz) < T_trail+1);
                ind_trail_max = max(inds_trail);
                N_bins_trail = ind_trail_max; 

                %Lump remaining bins
                Tu = T_liq_local;
                imax = length(inds_remaining); % imax = length(Tc);       
                imin = imax;
                
                remaining = true; %boolean - determines if there are possible cases to merge remaining
                while imax > ind_trail_max+1 

                    %Identify equantities to marge 
                    Tsz = Tcz_total(inds_remaining(1:imax), iz);


                    dT = Tu - Tsz; 
                    inds_merge = find(dT < dT_cutoff); %Indices to merge for this new bin
                    N_merge = max([length(inds_merge) - 2, 0]);
                    if N_merge > 0
                        imin = max([min(inds_merge), ind_trail_max]);
                    else
                        imin = imax;
                    end
                    

                    %Handle case of no mergers
                    if imin == Nt || N_merge == 0 %Move on to next bin and test
                        imax = imax - 1;
                        Tu = Tcz(inds_remaining(imax), iz);
                    else %Lump the bins 

                        %Define indices to merge
                        % if imax == Nt
                        %     imerged = (imin+1):(imax-1);
                        % else
                        imerged = (imin+1):(imax-1);  %exclude imax for merging
                        % end
                        imerged = inds_remaining(imerged);


                        %Mark the bin bounds to eliminate
                                % Tc_new(end+1,:) = Tcz(imin, iz:params.Nz);
                        Tcz(imerged, iz:params.Nz) = -1; %Set merged values to negative 1 to indicate 
    
                        %Log indexes correspondng to bin boundaries
                        

                        %Update bin inds
                        inds_bin(imerged) = -1;
    
    
                        %Set new Tu and imax
                        Tu = Tcz(inds_remaining(imin), iz); %Set the new upper limit
                        imax = imin;
                    end


                end
            
                %Determine inds remaining
                inds_remaining = find(Tcz(:,iz) > 0); 
                

                %Tc_new = [Tc_new; Tcz(inds_trail, iz)];
                N_bins(iz) = length(find(Tcz(:, iz) > -1));

                %Debug print
                if params.debug
                    fprintf('Merging Bins: im = %d; ix = %d; iz = %d; N_bins_orig = %d; N_bins_new = %d\n', im, ix, iz, Nt, N_bins(iz))
                end

                %Tcz_reduced{iz} = 1;
            end

            %Plot the final bin boundaries at each level
            if params.debug
                subplot(1,2,2);
                for iz = 1:params.Nz
                    z = params.mesh.volcell_cents(iz,2);
                    Ts = Tcz(:,iz); Ts = Ts(Ts>0);
                    zs = z .* ones(size(Ts));
                    plot(Ts, zs, 'r.', 'MarkerSize', 18); hold on;

                end
                close(debug_fig);

            end

            %Post-process results to determine which bins merge into which 
            N_bins = N_bins - 1;
            figure();
            
            bins_in = zeros(Nt-1, params.Nz);
            bins_in(:,1) = fliplr([1:(Nt-1)])';
            for iz = 1:(params.Nz) %Formerly + 1
                Tc = Tcz(:,iz);
                inds_merged = find(Tc(1:end-1) == -1);
                bins_merged = bins_in(inds_merged, iz);
                if iz ~= params.Nz
                    Tc_next = Tcz(:,iz+1);
                    ct = 1;
                    for ii = 1:(length(Tc_next)-1)
                        ind = length(Tc_next) - ii;
                        if Tc_next(ind) > -1
                            bins_in(ind,iz+1) = ct;
                            ct = ct + 1;
                        else
                            bins_in(ind,iz+1) = ct;
                        end   
                    end

                end

                %Debug plot
                
                if iz == 1
                    active_inds = 1:(Nt);
                    zs_plot = linspace(0, params.zms(iz), 20);
                elseif iz == params.Nz
                    active_inds = Tcz(:, iz-1) > 0;
                    zs_plot = linspace(params.zms(iz-1), params.reactor.H, 20);
                else
                    active_inds = Tcz(:, iz-1) > 0;
                    zs_plot = linspace(params.zms(iz-1), params.zms(iz), 20);
                end
                active_inds = find(active_inds);
                for i = 1:length(active_inds(1:end-1))
                    ia = active_inds(i);
                    Tz_plot = interp1(zs_all, Tzs{ia}{ix}, zs_plot);
                    subplot(1,2,1);
                    plot(Tz_plot, zs_plot, 'r-'); hold on
                    subplot(1,2,2);
                    plot(Tz_plot, zs_plot, 'r-'); hold on


                end
                

                subplot(1,2,1);
                grid on; grid minor; axis square;
                xlabel('Temperature (K)'); ylabel('z (m)')

                subplot(1,2,2);
                grid on; grid minor; axis square;
                xlabel('Temperature (K)'); ylabel('z (m)')
                xlim([1350, 1550]);
                
                Ts_z = 1;


            end

            N_cells(im) = sum(N_bins); %Number of cells required for the temperature
            


        end

        %Log results
        merge.Tcz = Tcz;
        merge.bins_in = bins_in;
        merge.bins_merged = bins_merged; 
        merged{im} = merge;


    end

    % %Create debug plot
    % if params.chars.debug
    %     for im = 1:Nm_max
    % 
    % 
    % 
    %     end
    % 
    % 
    % 
    % end



end