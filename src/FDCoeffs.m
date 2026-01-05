function mesh = FDCoeffs(mesh, deriv, opts)

%% Setup

    %Suppress warnings
    warning('off', 'all')

    %Run without inputs
    if nargin < 1
        PBM()
    end

    %Clean up
    %close all

    %Debug
    debug = true;

    %Create debug figure
    if debug
        debugfig = figure();
    end

    %Plot asestheics
    colors.focus = [1, 0.5, 0.5];
    colors.neighbors = [0.2, 0.2, 0.8];

%% Analysis

    %Pull values
    order = opts.order;
    direct_bias = opts.direct_bias;



if strcmp(mesh.type, '2DRegularCartesian')

    %Create discretization mat
    disc_mat = zeros(size(mesh.connmat)); %Matrix storing values 

    %Determine number of cells
    N_cells = mesh.N_cells;

    %
    xs_unique = unique(round(mesh.volcell_cents(:,1), 8));
    ys_unique = unique(round(mesh.volcell_cents(:,2), 8));
    N_xs = length(xs_unique);
    N_ys = length(ys_unique);

    %Check if order requested is greater than number of points in any one
    %direction
    if N_ys < order
        warning('Order exceeds mesh resolution. Automatically reducing.');
        order = min([N_xs, N_ys]) - 1;
    end
    
    F = zeros(1,order+1);
    F(deriv+1) = 1;

    %Iterate through cells
    for iv = 1:N_cells

        %Pull cell and neighbors
        volcell = mesh.volcells{iv};
        conns = mesh.connmat(iv,:);

        %Debug plot - cell and neighbors
        if debug
            cla
            plot(volcell.poly, 'FaceColor', colors.focus); hold on;
            for in = volcell.neigh_inds
                plot(mesh.volcells{in}.poly, 'FaceColor', colors.neighbors)
            end
        end

        %Find stencil cells
        [stencil_cells, stencil_dists] = FindStencilCells(iv, mesh, volcell, order, N_xs, N_ys, direct_bias);
        %stencil_cells = [xinds, yinds];
        stencil_poss = [mesh.volcell_cents(stencil_cells(:,1),1), mesh.volcell_cents(stencil_cells(:,2),2)];
        
        %Calculate coefficients for cell center
        cent = volcell.centroid;
        ycoeffs = CalcFDCoeffs(cent(2), stencil_poss(:,2), deriv);
        if mesh.Nxs == 1
            xcoeffs = zeros(size(ycoeffs));
        else
            xcoeffs = CalcFDCoeffs(cent(1), stencil_poss(:,1), deriv);
        end
        
        cent_FD_coeffs = [xcoeffs, ycoeffs];

        %Calculate coefficients for cell center - 2nd deriv
        ycoeffs = CalcFDCoeffs(cent(2), stencil_poss(:,2), 2);
        if mesh.Nxs == 1
            xcoeffs = zeros(size(ycoeffs));
        else
            xcoeffs = CalcFDCoeffs(cent(1), stencil_poss(:,1), 2);
        end
        
        cent_FD_coeffs2 = [xcoeffs, ycoeffs];

        %Calculate coefficients for each cell
        edge_FD_coeffs = cell(1,volcell.N_edges);
        for ie = 1:volcell.N_edges
            edge_cent = volcell.edge_cents(ie,:);

            %Calculate ycoeffs
            ycoeffs_D0 = CalcFDCoeffs(edge_cent(2), stencil_poss(:,2), 0);
            ycoeffs_D1 = CalcFDCoeffs(edge_cent(2), stencil_poss(:,2), deriv);

            %Calculate xcoeffs
            if mesh.Nxs == 1
                xcoeffs_D0 = zeros(size(ycoeffs));
                xcoeffs_D1 = zeros(size(ycoeffs));
            else
                xcoeffs_D0 = CalcFDCoeffs(edge_cent(1), stencil_poss(:,1), 0);
                xcoeffs_D1 = CalcFDCoeffs(edge_cent(1), stencil_poss(:,1), deriv);
            end
            

            %Log
            edge_FD_coeffs_D0{ie} = [xcoeffs_D0, ycoeffs_D0];
            edge_FD_coeffs_D1{ie} = [xcoeffs_D1, ycoeffs_D1];

        end

        %Log results
        mesh.volcells{iv}.stencil_cells = stencil_cells;
        mesh.volcells{iv}.stencil_dists = stencil_dists;
        mesh.volcells{iv}.cent_FD_coeffs_D1 = cent_FD_coeffs;
        mesh.volcells{iv}.cent_FD_coeffs_D2 = cent_FD_coeffs2;
        mesh.volcells{iv}.edge_FD_coeffs_D1 = edge_FD_coeffs_D1;
        mesh.volcells{iv}.edge_FD_coeffs_D0 = edge_FD_coeffs_D0;

    end



elseif strcmp(mesh.type, 'Arbitrary')
    error('Mesh Type not Implemented')

elseif strcmp(mesh.type, '1D')

    %Assumes points are ordered smallest to largest

    %Determine number of volumes
    N_cells = length(mesh.xsc);

    %Calculate cell centered values
    for ic = 1:length(mesh.xsc)

        %Allocate stencil vector
        stencil_inds = zeros(1,opts.order+1);
        stencil_inds(1) = ic; %Always use current cell 

        %Calculate nearest neighbors
        xc = mesh.xsc(ic);
        xb = mesh.xsb(ic);
        xsc_dists = mesh.xsc - xc;
        belows = xsc_dists < 0;
        N_below = length(find(belows));
        aboves = xsc_dists > 0;
        N_above = length(find(aboves));
        inds_below = 1:ic-1;
        inds_above = ic+1:N_cells;
        

        if opts.direct_bias == 1 %Forward differences

            %Use as many larger values as possible up to order requested
            Ni_above = min([N_above, order]);
            Ni_below = order - (Ni_above);

            %Define above inds
            if Ni_above > 0
                stencil_inds(2:(1+Ni_above)) = inds_above(1:Ni_above);
            end

            %Define below inds
            if Ni_below > 0
                stencil_inds(2+Ni_above:end) = inds_below((N_below - Ni_below + 1):end);
            end
           
        elseif opts.direct_bias == -1 %Backwards differences

            %Use as many larger values as possible up to order requested
            Ni_below = min([N_below, order]);
            Ni_above = order - (Ni_below);

            %Define above inds
            if Ni_above > 0
                stencil_inds(2:(1+Ni_above)) = inds_above(1:Ni_above);
            end

            %Define below inds
            if Ni_below > 0
                stencil_inds(2+Ni_above:end) = inds_below((N_below - Ni_below + 1):end);
            end

        else %Centered differences - try for at least one point on each side

            %Find nearest points
            Ni_below_ideal = floor((order+1)/2); %Bias to below because of log scaling
            Ni_above_ideal = floor((order)/2);

            %Account for edges
            if N_below < Ni_below_ideal
                Ni_below = N_below;
                Ni_above = order - Ni_below;
            elseif N_above < Ni_above_ideal
                Ni_above = N_above;
                Ni_below = order - Ni_above;
            else
                Ni_below = Ni_below_ideal;
                Ni_above = Ni_above_ideal;
            end
            
            %Define above inds
            if Ni_above > 0
                stencil_inds(2:(1+Ni_above)) = inds_above(1:Ni_above);
            end

            %Define below inds
            if Ni_below > 0
                stencil_inds(2+Ni_above:end) = inds_below((N_below - Ni_below + 1):end);
            end

            
        end

        %Reorder stencil cells
        stencil_inds = sort(stencil_inds);
        stencil_xcs = mesh.xsc(stencil_inds);
        coeffs = CalcFDCoeffs(xc, stencil_xcs, deriv);

        %Log results
        stencil_inds_c(:,ic) = stencil_inds;
        stencil_xcs_c(:,ic) = stencil_xcs;
        coeffs_c(:,ic) = coeffs;

        %Caclulate bracket coefficients
        coeffs_b(:,ic) = CalcFDCoeffs(xb, stencil_xcs, deriv);
        if ic == N_cells %Do final bracket position
            coeffs_b(:,ic+1) = CalcFDCoeffs(mesh.xsb(ic+1), stencil_xcs, deriv);
        end

    end

    

    %Calculate interface values (brackets)
    mesh.stencil_inds_c = stencil_inds_c;
    mesh.stencil_xcs_c = stencil_xcs_c;
    mesh.coeffs_c = coeffs_c;
    mesh.coeffs_b = coeffs_b;


end













end


function ClassifyNeighborDir()


end


function [stencil_cells, stencil_dists] = FindStencilCells(iv, mesh, volcell, order, N_xs, N_ys, direct_bias)

        %Specify cells to search 
        other_inds = [1:(iv-1), (iv+1):mesh.N_cells];

        %Calculate distances
        dist_vecs = abs(mesh.volcell_cents - volcell.centroid);
        [y_dists, y_dist_sort_inds] = sort(dist_vecs(:,2)); 
        [x_dists, x_dist_sort_inds] = sort(dist_vecs(:,1));

        %Determine nearest cells in the x direction
        if N_xs > 1
            xs_nearest_inds = y_dist_sort_inds(1:N_xs); %Take nearest points
            xs_nearest = mesh.volcell_cents(xs_nearest_inds, 1); 
            xs_dists = abs(xs_nearest - volcell.centroid(1));
            [xs_dists_sorted, xs_dists_sorted_inds] = sort(xs_dists); %Sort
            xs_dists_sorted_inds = xs_nearest_inds(xs_dists_sorted_inds);
            xinds = xs_dists_sorted_inds(1:(order+1));
            xposs = mesh.volcell_cents(xinds,1);
            [xposs, align_inds] = sort(xposs);
            xdists = xposs - volcell.centroid(1);
            xinds = xinds(align_inds); % first = most left, last = most right
        else
            xdists = zeros(order+1,1);
            xinds = iv * ones(order+1,1);
        end

        %Determine nearest cells in the y direction
        ys_nearest_inds = x_dist_sort_inds(1:N_ys);
        ys_nearest = mesh.volcell_cents(ys_nearest_inds,2);
        ys_dists = abs(ys_nearest - volcell.centroid(2));
        N_neg = length(find(ys_dists < 0));
        N_pos = length(find(ys_dists > 0));
        % if direct_bias == 0 || N_neg < (order+1) %Will default to centered differences
        [ys_dists_sorted, ys_dists_sorted_inds] = sort(ys_dists); %Sort
        ys_dists_sorted_inds = ys_nearest_inds(ys_dists_sorted_inds);
        yinds = ys_dists_sorted_inds(1:(order+1));
        % elseif direct_bias == -1 %Will use backwards  differences wherever possible   -MODIFY FOR EDGE CALCULATIONS TO STILL USE CELL CENTER
        %     neg_inds = ys_dists <= 0
        % 
        % 
        % 
        % end


        yposs = mesh.volcell_cents(yinds,2);
        [yposs, align_inds] = sort(yposs);
        ydists = yposs - volcell.centroid(2);
        yinds = yinds(align_inds);

        %Define outputs
        stencil_cells = [xinds, yinds];
        stencil_dists = [xdists, ydists];


end

function coeffs = CalcFDCoeffs(x0, xos, deriv)
    N_pts = length(xos);
    A = zeros(N_pts);
    for i = 1:length(xos)
        for j = 1:length(xos)
            A(i,j) = ((xos(j) - x0).^(i-1))./factorial(i-1);
        end
    end

    %Solve system
    F = zeros(N_pts,1);
    F(deriv+1) = 1;
    coeffs = A\F;

end