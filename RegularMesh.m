function mesh = RegularMesh(xx, yy, zz, options)

%% Setup

    %Code time
    t_in = cputime;

    %Enable run without changing files
    if nargin < 1
        PBM();
    end

    shape = options.shape;

    %Clean up
    %close all

    %Debug
    debug = false;

    if nargin < 3
        N_dims = 2;
        zz = zeros(size(xx));
        [Ny, Nx] = size(xx);
        Nz = 1;
    else
        N_dims = 3;
        [Ny, Nx, Nz] = size(xx);
    end

%% Mesh Generation
%
%    Assume that the points correspond to rectangular array

    %Store input points
    mesh.xx = xx;
    mesh.yy = yy;
    mesh.zz = zz;


    %Create cells array
    N_cells = (Nx-1) .* (Ny-1);
    volcells = cell(N_cells,1);

    %Define cells
    it = 1;
    for ix = 1:(Nx-1)
        for iy = 1:(Ny-1)
            
            %Define cell index
            volcell.ind = it;

            %Define cell points
            volcell.pts = [xx(iy, ix), yy(iy, ix);
                        xx(iy, ix+1), yy(iy, ix+1);
                        xx(iy+1, ix+1), yy(iy+1,ix+1);
                        xx(iy+1, ix), yy(iy+1, ix)];
            volcell.xs = volcell.pts(:,1); volcell.ys = volcell.pts(:,2);
            volcell.N_pts = length(volcell.xs); 

            
          
            %Define cell edges
            volcell.edges = [1, 2;
                          2, 3;
                          3, 4;
                          4, 1];
    
            %Log number of edges
            [volcell.N_edges, ~] = size(volcell.edges);

            %Define cell polyshape
            volcell.poly = polyshape(volcell.xs, volcell.ys);

            %Define cell centroid
            [cx, cy] = centroid(volcell.poly);
            volcell.centroid = [cx, cy];

            %Calculate cell area and volume
            volcell.area = area(volcell.poly);
            if strcmp(shape, 'Rect_2D')
                volcell.vol = volcell.area; %assume 1 m depth
            elseif strcmp(shape, 'Cylinder') %Assume 1st dimension is radius
                [volcell.vol, volcell.surfareas] = CalcCylindricalVolume(volcell);
            end

            %Define edge centers
            volcell.edge_cents = CalcEdgeCenters(volcell.pts, volcell.edges);

            %Calculate normal vector for each edge
            volcell.edge_ns = CalcNormVectors(volcell);

            %Allocate neigh_inds
            volcell.neigh_inds = [];

            %Determine orientation of edges
            volcell.edge_orients = [];

            %Store cell centroids
            mesh.volcell_cents(it,:) = volcell.centroid;


            %Log cell and update counter
            volcells{it} = volcell;
            it = it + 1;

        end
    end

    %Store cells
    mesh.volcells = volcells;

    %Check how many 
    mesh.Nxs = length(unique(mesh.volcell_cents(:,1)));
    mesh.Nys = length(unique(mesh.volcell_cents(:,2)));

    %Debug plot
    if debug & false
        PlotMesh(mesh)
    end

    %Identify neighbors
    connmat = zeros(N_cells, N_cells);
    for ic = 1:N_cells
        cell1 = mesh.volcells{ic};
        ind = cell1.ind;
        connmat(ic,ic) = 1;
        %edges = mesh.volcells{ic}.edges;
        %N_edges = mesh.volcells{ic}.N_edges;
        neigh_inds = [];
        inds_other = [1:(ind-1), (ind+1):N_cells];
        
        %Iterate through all other cells
        for io = inds_other
            cell2 = mesh.volcells{io};

            %Check if other cell has already determined this cell is a neighbor
            if any(cell2.neigh_inds == ic)
                connmat(ic,io) = 1;
                neigh_inds(end+1) = io;
            elseif io > ic %only test cells furhter down line
                neighbors = CompareCells(cell1, cell2);
                if neighbors
                    neigh_inds(end+1) = io;
                    connmat(ic,io) = 1;
                end
            end
        end




        % %Iterate through each edge
        % for ie = 1:N_edges
        %     pt1 = mesh.volcells{ic}.pts(edges(ie, 1), :);
        %     pt2 = mesh.volcells{ic}.pts(edges(ie, 2), :);
        %     edge_mat = [pt1; pt2];
        % 
        % 
        %     %Iterate through all other cells
        %     for io = inds_other
        % 
        % 
        % 
        % 
        % 
        % 
        %         % if connmat(ic, io) == 1 && ~any(neigh_inds == io) %already connected, just add to storage vec
        %         %     neigh_inds(end+1) = io;
        %         % else %Not already shown to be connected 
        %         % 
        %         %     edges_o = mesh.volcells{io}.edges;
        %         %     N_edges_o = mesh.volcells{io}.N_edges;
        %         % 
        %         %     %Check all edges
        %         %     for ie2 = 1:N_edges_o
        %         %         pt1 = mesh.volcells{io}.pts(edges(ie2,1), :);
        %         %         pt2 = mesh.volcells{io}.pts(edges(ie2,2), :);
        %         %         edge_mat_o = [pt1; pt2];
        %         % 
        %         %         %Check if edge is the same 
        %         %         equiv = CompareEdges(edge_mat, edge_mat_o);
        %         % 
        %         %         %Log neighborsd
        %         %         if equiv
        %         %             neigh_inds(end+1) = io;
        %         %             connmat(ic,io) = 1;
        %         %         end
        %         % 
        %         %     end
        %         % 
        %         % end
        %     end
        % 
        % 
        % 
        % 
        % 
        % end

        %Log neighbor inds
        mesh.volcells{ic}.neigh_inds = neigh_inds;
        
    end
    
    %Calculate normal vectors


    %Calulate height of each cell
    mesh.hs = mesh.yy(2:end,1) - mesh.yy(1:end-1,1);

    %Define outputs
    mesh.connmat = connmat;
    mesh.type = '2DRegularCartesian';
    mesh.N_cells = length(mesh.volcells);

    %Calculate and report cpu time
    t_out = cputime;
    t_run = abs(t_out - t_in);
    fprintf('Mesh Generation Complete (t = %0.2f s)\n', t_run);


end


function centers = CalcEdgeCenters(pts_in, conns)

    [N_pts, ~] = size(pts_in);
    centers = zeros(N_pts, 2);
    for i = 1:N_pts
        edge = conns(i,:);
        pts = pts_in(edge',:);
        centers(i,:) = sum(pts)./2;
    

    end
end

function [vol, surfareas] = CalcCylindricalVolume(volcell)

    %Pull important values
    pts = volcell.pts;
    N_edges = volcell.N_edges;

    %Iterate through edges and represent as a line
    lines = cell(volcell.N_pts,1);
    ms = zeros(size(lines));
    bs = ms;
    for i = 1:volcell.N_pts
        edge = volcell.edges(i,:);
        edge_pts = volcell.pts(edge, :);
        if edge_pts(1,1) == edge_pts(2,1)
            ms(i) = NaN;
            bs(i) = NaN;
            lines{i} = @(x) NaN;
        else
            m = (edge_pts(2,2) - edge_pts(1,2))./(edge_pts(2,1) - edge_pts(1,1));
            b = edge_pts(1,2) - m .* edge_pts(1,1);
            lines{i} = @(x) m.*x + b;
            ms(i) = m;
            bs(i) = b;
        end
    end

    %Evaluate integral - use first edge as baseline 
    intsum = 0;

    %Area of left side (-r) triangle
    if pts(4,1) < pts(1,1) && ~isnan(lines{4}(0))%Integrate l3 - l2

        %Volume
        xi = pts(4,1);
        xf = pts(1,1);
        Hfunc = @(x) lines{3}(x) - lines{4}(x);
        Vfunc = @(x) 2 .* pi .* x * Hfunc(x);
        intsum = intsum + integral(Vfunc, xi, xf);

        %Area

    else


    end
    
    if pts(3,1) > pts(2,1) && ~isnan(lines{2}(0))
        xi = pts(2,1);
        xf = pts(3,1);
        Hfunc = @(x) lines{3}(x) - lines{2}(x);
        Vfunc = @(x) 2 .* pi .* x .* Hfunc(x);
        intsum = intsum + integral(Vfunc, xi, xf);
    end

    %Area of main segment
    xi = pts(1,1);
    xf = pts(2,1);
    Hfunc = @(x) lines{3}(x) - lines{1}(x);
    Vfunc = @(x) 2 .* pi .* x .* Hfunc(x);
    intsum = intsum + integral(Vfunc, xi, xf);

    %Calculate Edge surface areas 
    mcorrect = sqrt(ms(1)^2 + 1);
    A1func = @(x) 2*pi*mcorrect*x;
    areas(1) = integral(A1func, pts(1,1), pts(2,1));

    mcorect = sqrt(ms(3)^2 + 1);
    A3func = @(x) 2*pi*mcorrect*x;
    areas(3) = integral(A3func, pts(4,1), pts(3,1));

    if pts(3,1) > pts(2,1) && ~isnan(lines{2}(0))
        mcorrect = sqrt(ms(2)^2 + 1);
        A2func = @(x) 2*pi*mcorrect*x;
        areas(2) = integral(A2func, pts(2,1), pts(3,1));
    else
        areas(2) = (pts(3,2) - pts(2,2)) * 2*pi*pts(2,1); %Vertical case
    end

    if pts(4,1) < pts(1,1) && ~isnan(lines{4}(0))
        mcorrect = sqrt(ms(4)^2 + 1);
        A4func = @(x) 2*pi*mcorrect*x;
        areas(4) = integral(A2func, pts(4,1), pts(1,1));
    else
        areas(4) = (pts(4,2) - pts(1,2)) * 2*pi*pts(1,1); %Vertical case
    end

    %Assign output
    vol = intsum;
    surfareas = areas;


end


function equiv = CompareEdges(edge1, edge2)

    
    edge_match = true;
    for i1 = 1:2
        pt1 = edge1(i1,:);
        for i2 = 1:2
            pt2 = edge2(i2,:);

            %Check if points match
            
            pts_match = all(pt1 == pt2);

            %Exit if match
            if pts_match
                break
            end

        end

        edge_match = edge_match && pts_match;

    end



    equiv = edge_match;
end


function neighbors = CompareCells(cell1, cell2)

    %Define output default
    neighbors = false;

    %Pull points
    pts1 = cell1.pts;
    pts2 = cell2.pts;

    N_pts1 = cell1.N_pts;
    N_pts2 = cell2.N_pts;

    %Check if two points are shared
    pts_matching = 0;
    for i1 = 1:N_pts1

        %Pull point from first cell
        pt1 = pts1(i1,:);

        for i2 = 1:N_pts2

            %Pull point from second cell
            pt2 = pts2(i2,:);

            %Check if points match
            if pt1 == pt2
                pts_matching = pts_matching + 1;
                if pts_matching >= 2
                    neighbors = true;
                    break
                end
            end

        end
    end
end

function edge_ns = CalcNormVectors(volcell)

    %Allocate 
    edge_ns= zeros(volcell.N_edges, 2);

    %Iterate through edges
    for ie = 1:volcell.N_edges
        inds = volcell.edges(ie,:);
        pt1 = volcell.pts(inds(1), :);
        pt2 = volcell.pts(inds(2), :);


        edge_vec = pt2 - pt1;
        edge_len = norm(edge_vec);
        edge_ns(ie,:) = [edge_vec(2), -edge_vec(1)]/edge_len;

    end

end