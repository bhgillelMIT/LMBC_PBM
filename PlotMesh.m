function PlotMesh(mesh, params)

    %Plot settings
    ms = 12;
    fs = 18;
    lw = 2;

    %Optional inputs
    if nargin < 2
        cellcolor = [0.7, 0.7, 0.7];
        plotcents = true;
        plotedgecents = true;
    else
        cellcolor = params.cellcolor;
        plotcents = params.plotcents;
        plotedgecents = params.plotedgecents;
    end

    %Create figure
    figure();

    %Plot cells
    N_cells = length(mesh.volcells);
    for ic = 1:N_cells

        %Plot cells
        volcell = mesh.volcells{ic};
        plot(volcell.poly, 'FaceColor', cellcolor); hold on;

        %Plot centroids
        if plotcents
            plot(volcell.centroid(1), volcell.centroid(2), 'r.', 'MarkerSize', ms)
        end

        %Plot edge centers
        if plotedgecents
            [N_edges, ~] = size(volcell.edges);
            for ie = 1:N_edges
                plot(volcell.edge_cents(ie, 1), volcell.edge_cents(ie, 2), 'b.', 'MarkerSize', ms-4);
            end
        end
    end

    title('Mesh');
    axis equal

    


end