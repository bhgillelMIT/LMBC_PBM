function [y0, params, N_volumes_total] = InitializeVolumes(params, mesh, disc, Fb_i)

%% Setupp

    if nargin < 1
        PBM_v3()
    end

    %Settings
    dT_cutoff = 10; %K
    dX_cutoff = 0.02; % percent - 


    %Plot aesthetis
    ms = 18;



%% For bubbles that converge to equilibrium temperature and conversion 

    % %Create output cell array
    % merged = cell(1,params.Nms);
    % 
    % %Iterate through mass
    % 
    % 
    % for im = 1:disc.Nms
    % 
    %     %Pull characteristics
    %     Tzs = params.T_z.chars.Ts_zs{im};
    %     Tcs = params.T_z.chars.Ts_cs{im};
    % 
    %     %Determine number
    %     Nx = length(Tcs{1});
    % 
    %     %Define boundaries
    %     %bins_red = params.T_liq:-dT_cutoff: 
    % 
    % 
    %     % %Debug plot
    %     % if params.debug
    %     %     for ix = 1:Nx
    %     % 
    %     %     end
    %     % 
    %     % 
    %     % end
    % 
    % 
    % 
    %     %Iterate through initial conversion cases
    %     for ix = 1:1 %1:Nx
    % 
    %         %Debug plot
    %         if params.debug
    %             debug_fig = figure();
    %             subplot(1,2,1);
    %             for it = 1:length(Tcs)
    %                 plot(Tzs{it}{ix}, linspace(0, params.reactor.H, length(Tzs{it}{ix})), 'r-'); hold on;
    %                 plot(Tcs{it}{ix}, params.mesh.volcell_cents(:,2), 'r.', 'MarkerSize', ms); hold on;
    %             end
    %             title_str = sprintf('Nm = %d; Nx = %d;', im, ix);
    %             title(title_str);
    %         end
    % 
    % 
    %         %Pull current bin temperatures
    %         Nt = length(Tcs);
    %         Tcz = zeros(Nt, params.Nz);
    %         for it = 1:Nt
    %             Tcz(it,:) = Tcs{it}{ix}';
    %         end
    % 
    % 
    %         %Append liquid temperature to Tcz
    %         for iz = 1:params.Nz
    %             Tcz(Nt+1, iz) = params.T_Lz(params.mesh.volcell_cents(iz, 2));
    %         end
    %         Nt = Nt + 1;
    %         Tcz_total = Tcz; %Copy of the original bins
    % 
    %         %Allocate storage for new bounds
    %         inds_bin = 1:Nt;
    %         N_bins = zeros(1, params.Nz);
    % 
    %         %Iterate through Z-levels 
    %         inds_remaining = inds_bin;
    %         for iz = 1:params.Nz
    % 
    %             %Isolate this levels temperature bounds and count bins
    %             Tc = Tcz(:,iz);
    %             %Nt = Nt + 1;
    %             N_bins_i = length(Tcz(:,iz))-1;
    % 
    %             %Determine initial bounds
    %             T_trail = Tcz(1,iz);
    %             T_liq_local = params.Tsz(iz);
    %             Tc_new = T_liq_local.*ones(1, params.Nz);
    % 
    %             %Determine how bins are trailing points
    %             inds_trail = find(Tcz_total(:,iz) < T_trail+1);
    %             ind_trail_max = max(inds_trail);
    %             N_bins_trail = ind_trail_max; 
    % 
    %             %Lump remaining bins
    %             Tu = T_liq_local;
    %             imax = length(inds_remaining); % imax = length(Tc);       
    %             imin = imax;
    % 
    %             remaining = true; %boolean - determines if there are possible cases to merge remaining
    %             while imax > ind_trail_max+1 
    % 
    %                 %Identify equantities to marge 
    %                 Tsz = Tcz_total(inds_remaining(1:imax), iz);
    % 
    % 
    %                 dT = Tu - Tsz; 
    %                 inds_merge = find(dT < dT_cutoff); %Indices to merge for this new bin
    %                 N_merge = max([length(inds_merge) - 2, 0]);
    %                 if N_merge > 0
    %                     imin = max([min(inds_merge), ind_trail_max]);
    %                 else
    %                     imin = imax;
    %                 end
    % 
    % 
    %                 %Handle case of no mergers
    %                 if imin == Nt || N_merge == 0 %Move on to next bin and test
    %                     imax = imax - 1;
    %                     Tu = Tcz(inds_remaining(imax), iz);
    %                 else %Lump the bins 
    % 
    %                     %Define indices to merge
    %                     % if imax == Nt
    %                     %     imerged = (imin+1):(imax-1);
    %                     % else
    %                     imerged = (imin+1):(imax-1);  %exclude imax for merging
    %                     % end
    %                     imerged = inds_remaining(imerged);
    % 
    % 
    %                     %Mark the bin bounds to eliminate
    %                             % Tc_new(end+1,:) = Tcz(imin, iz:params.Nz);
    %                     Tcz(imerged, iz:params.Nz) = -1; %Set merged values to negative 1 to indicate 
    % 
    %                     %Log indexes correspondng to bin boundaries
    % 
    % 
    %                     %Update bin inds
    %                     inds_bin(imerged) = -1;
    % 
    % 
    %                     %Set new Tu and imax
    %                     Tu = Tcz(inds_remaining(imin), iz); %Set the new upper limit
    %                     imax = imin;
    %                 end
    % 
    % 
    %             end
    % 
    %             %Determine inds remaining
    %             inds_remaining = find(Tcz(:,iz) > 0); 
    % 
    % 
    %             %Tc_new = [Tc_new; Tcz(inds_trail, iz)];
    %             N_bins(iz) = length(find(Tcz(:, iz) > -1));
    % 
    %             %Debug print
    %             if params.debug
    %                 fprintf('Merging Bins: im = %d; ix = %d; iz = %d; N_bins_orig = %d; N_bins_new = %d\n', im, ix, iz, Nt, N_bins(iz))
    %             end
    % 
    %             %Tcz_reduced{iz} = 1;
    %         end
    % 
    %         %Plot the final bin boundaries at each level
    %         if params.debug
    %             subplot(1,2,2);
    %             for iz = 1:params.Nz
    %                 z = params.mesh.volcell_cents(iz,2);
    %                 Ts = Tcz(:,iz); Ts = Ts(Ts>0);
    %                 zs = z .* ones(size(Ts));
    %                 plot(Ts, zs, 'r.', 'MarkerSize', 18); hold on;
    % 
    %             end
    %             close(debug_fig);
    % 
    %         end
    % 
    %         %Post-process results to determine which bins merge into which 
    %         N_bins = N_bins - 1;
    %         bins_in = zeros(Nt-1, params.Nz);
    %         bins_in(:,1) = fliplr([1:(Nt-1)])';
    %         for iz = 1:params.Nz
    %             Tc = Tcz(:,iz);
    %             inds_merged = find(Tc(1:end-1) == -1);
    %             bins_merged = bins_in(inds_merged, iz);
    %             if iz ~= params.Nz
    %                 Tc_next = Tcz(:,iz+1);
    %                 ct = 1;
    %                 for ii = 1:(length(Tc_next)-1)
    %                     ind = length(Tc_next) - ii;
    %                     if Tc_next(ind) > -1
    %                         bins_in(ind,iz+1) = ct;
    %                         ct = ct + 1;
    %                     else
    %                         bins_in(ind,iz+1) = ct;
    %                     end   
    %                 end
    % 
    %             end
    % 
    %         end
    % 
    % 
    % 
    %     end
    % 
    %     %Log results
    %     merge.Tcz = Tcz;
    %     merge.bins_in = bins_in;
    %     merge.bins_merged = bins_merged; 
    %     merged{im} = merge;
    % 
    % 
    % end




%% Determine the number of volumes


    %Pull merged 
    if params.heat.active
        chars = params.T_z.chars;
    end


    N_volumes_total = 0;
    N_volumes_layer = zeros(disc.Nz, 1);
    
    cellinds = [];
    xinds = [];
    zinds = [];
    Tinds = [];
    Xinds = [];


    %Iterate through each z layer
    for iz = 1:disc.Nz

        if params.heat.active & params.react.active
            
 
                for im = 1:disc.Nms

                    


                    NTs = length(params.T_z.chars.Ts_zs{im});
                    NXs = length(params.T_z.chars.Ts_zs{im}{1});
                    Npts = NTs*NXs;
                    N_volumes_layer(iz) = N_volumes_layer + Npts;
                end
         
        elseif params.heat.active

            %Allocate
            cellinds_layer = []; %[1:mesh.N_cells], disc.Nms, 1); params.cellinds = params.cellinds(:); %Specifies which spatial cell the volume maps to
            xinds_layer = []; %repmat([1:disc.Nms]', mesh.N_cells, 1);
            zinds_layer = []; %repmat([1:disc.Nz], disc.Nms, 1); params.zinds = params.zinds(:);
            Tinds_layer = [];
            Xinds_layer = [];

            %Iterate through each representative size
            for im = 1:disc.Nms

                


                %Determine how many temperatures are relevant
                merged = chars.merged{im};
                bins_in = merged.bins_in(:,iz);
                NTs = max(bins_in);
                NXs = 1;
                Npts = NTs*NXs;
                N_volumes_layer(iz) = N_volumes_layer(iz) + Npts;

                %Define indices
                cellinds_size = iz .* ones(Npts, 1); %[1:mesh.N_cells], disc.Nms, 1); params.cellinds = params.cellinds(:); %Specifies which spatial cell the volume maps to
                xinds_size = im .* ones(Npts, 1); %repmat([1:disc.Nms]', mesh.N_cells, 1);
                zinds_size = iz .* ones(Npts, 1); %repmat([1:disc.Nz], disc.Nms, 1); params.zinds = params.zinds(:);
                Tinds_size = 1:NTs; Tinds_size = Tinds_size';
                Xinds_size = ones(Npts, 1);

                cellinds_layer = [cellinds_layer; cellinds_size];
                xinds_layer = [xinds_layer; xinds_size];
                zinds_layer = [zinds_layer; zinds_size];
                Tinds_layer = [Tinds_layer; Tinds_size];
                Xinds_layer = [Xinds_layer; Xinds_size];


            end

            % cellinds_layer = [cellinds_layer; cellinds_size];
            % xinds_layer = [xinds_layer; xinds_size];
            % zinds_layer = [zinds_layer; zinds_size];
            % Tinds_layer = [Tinds_layer; Tinds_size];
            % Xinds_layer = [Xinds_layer; Xinds_size];
        


        elseif params.react.active 
                x = 1

        else %No heat or reaction - default case 

            %Define number of volumes and indices
            N_volumes_layer(iz) = disc.Nms;
            cellinds_layer = iz .* ones(disc.Nms, 1); %[1:mesh.N_cells], disc.Nms, 1); params.cellinds = params.cellinds(:); %Specifies which spatial cell the volume maps to
            xinds_layer = 1:1:disc.Nms; xinds_layer = xinds_layer'; %repmat([1:disc.Nms]', mesh.N_cells, 1);
            zinds_layer = iz .* ones(disc.Nms, 1); %repmat([1:disc.Nz], disc.Nms, 1); params.zinds = params.zinds(:);
            Tinds_layer = ones(disc.Nms, 1);
            Xinds_layer = ones(disc.Nms, 1);
            
            %bottom_inds = params.cents_y < min(mesh.volcell_cents(:,2)) + 1E-6;



        end
    
        %Log layer values
        cellinds = [cellinds; cellinds_layer];
        xinds = [xinds; xinds_layer];
        zinds = [zinds; zinds_layer];
        Tinds = [Tinds; Tinds_layer];
        Xinds = [Xinds; Xinds_layer];



    end

    %Identify bottom indices
    bottom_inds = zinds == 1;


    %Store indices
    params.cellinds = cellinds;
    params.xinds = xinds;
    params.zinds = zinds;
    params.Tinds = Tinds;
    params.Xinds = Xinds;




    %Calculate total number of volumes
    N_volumes_total = sum(N_volumes_layer);
    N_volumes_total_nomerge = N_volumes_layer(1) * params.Nz;
    



    %Calculate total number of volumes

%% Allocate those volumes

    %Layout of vector
    % - 4 levels
    %    - level 1 - z layer 
    %    - level 2 - mass (representative

    %Create initial "volumes" and specify boundary conditions
    y0 = zeros(N_volumes_total, 1);
    bottom_inds = params.cents_y < min(mesh.volcell_cents(:,2)) + 1E-6;
    if params.sol.orifice_BC_type == 1
        Fi_bottom = Fb_i; %repmat(Fb_i, 1,);
        y0(bottom_inds) = Fi_bottom;
    end

    %Create storage vector 
    %y0 = zeros(N_volumes_total, 1);

    %Specify initial values
    
 
    %Define reference cell values
    Vcells_rep = zeros(N_volumes_total,1);
    mms_rep = zeros(N_volumes_total,1);

    for iz = 1:params.Nz

        %Store cell volume (spatial)
        if iz == 1
            Vcells_rep(1:N_volumes_layer(1)) = params.Vcells(iz);
        else
            Vcells_rep(sum(N_volumes_layer(1:iz-1)):sum(N_volumes_layer(1:iz))) = params.Vcells(iz);
        end

        

    end

    %Iterate 
    for im = 1:params.Nms
        mms_rep(xinds == im) = params.mms(im);
    end

    params.Vcells_rep = Vcells_rep; %repmat(params.Vcells, 1, params.Nms); params.Vcells_rep = params.Vcells_rep';
    %params.Vcells_rep = params.Vcells_rep(:);
    params.mms_rep = mms_rep; %repmat(params.mms', params.Nz, 1);



end