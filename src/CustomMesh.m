function mesh = CustomMesh(xx, yy, zz, shape)

%% Setup

    %
    if nargin < 1
        PBM();
    end

    %Clean up
    close all

    %Debug
    debug = true;

    if nargin < 3
        N_dims = 2;
        zz = zeros(size(xx));
        [Nx, Ny] = size(xx);
        Nz = 1;
    else
        N_dims = 3;
        [Nx, Ny, Nz] = size(xx);
    end

%% Identify closest neighbors

    %Debug plot
    if debug
        scatter3(xx, yy, zz)
    end

    %Identify range of each domain
    xbs = [min(min(xx)), max(max(xx))];
    ybs = [min(min(yy)), max(max(yy))];
    zbs = [min(min(zz)), max(max(zz))];

    %Identify unique elements
    xs_unique = unique(xx);
    ys_unique = unique(yy);
    zs_unique = unique(zz);

    %Print update
    fprintf('Generating connectivity matrix');

    %Generate index mat
    N_pts = numel(xx);
    ind_vec = [1:N_pts];
    ind_mat = reshape(ind_vec, Nx, Ny, Nz);

    %Iterate throguh points
    
    it = 1;
    


    for ix = 1:Nx
        for iy = 1:Ny
            for iz = 1:Nz
        
                

                %Pull values
                x = xx(ix, iy, iz);
                y = yy(ix, iy, iz);
                z = zz(ix, iy, iz);
        
                
        
                %Find closest neighbors
                if N_dims == 2

                    %Calculate distances
                    dists_x = sqrt((xx(:,:,:) - x).^2);
                    dists_y = sqrt((yy(:,:,:) - y).^2);
                    dists = sqrt( (xx - x).^2 + (yy - y).^2 + (zz - z).^2);

                    %Determine if on boundaries
                    on_xb = any(x == xbs); %on x boundaries
                    on_yb = any(y == ybs); %on y boundaries
                    on_bb = on_xb && on_yb; %on both boundaries 

                    %Sort distances
                    dists_x_sorted = sort(dists_x(:));
                    dists_y_sorted = sort(dists_y(:));
                    [dists_sorted, inds_sorted] = sort(dists(:));

                    %
                    if on_bb %corner - find one x and one y
                        
                        closest_pts = 1;

                        closest_x = min(abs(dists));
                        closest_y = 1;

                        %Find nearest x in pos dir - min nonzero xdist and min
                        %ydist

                        %Find nearest x in neg dir

                        %Find nearest y in pos dir



                        %Find nearest y in neg dir


                        xns = 1;
                        yns = 1;
                    elseif on_xb
                        xns = 1;
                        yns = 1;
                    elseif on_yb
                        xns = 1;
                        yns = 1;
                    else %interior
                        xns = 1;
                        yns = 1;
                    end

                    

                    %Isolate four nearest points
                    dists_sorted = sort(dists(:));
                    dists_nearest = dists_sorted(1:5);
                    near_pts = 1;
        
                    %Find nearest neighbors in x
                    if any(x == xbs) %edge
                        ins = 1;
                        xns = 1;
                    else %not edge
        
                    end
        
        
                    %Find nearest neighbors in y
                else
                    dists_z = sqrt((zz(:,:,iz) - z).^2);
        
        
                end

                %Iteration counter
                it = it + 1;

            end
        end


    end

    %Define output
    mesh = []

end

function FindInds()


end