function mat = PBM_Mat(mesh, opts)



    %Pull options
    z_dir = opts.z_dir; 

    %Allocate Mat
    mat = zeros(mesh.N_cells);

    %Iterate through cells
    for iv = 1:mesh.N_cells

        %Pull cell and relevant values
        volcell = mesh.volcells{iv};
        cent = volcell.centroid;
        ystencil = volcell.stencil_cells(:,2);
        cent_coeffs_D1 = round(volcell.cent_FD_coeffs_D1(:,2), 8);
        edge_coeffs_D0 = volcell.edge_FD_coeffs_D0;
        edge_coeffs_D1 = volcell.edge_FD_coeffs_D1;

        %F Derivative - Using cell centered FD
        for i = 1:length(ystencil)
            io = ystencil(i);
            val = cent_coeffs_D1(i);
            mat(iv,io) = val;
            x=1;

            

        end
        z_edges = 1;

    end

end