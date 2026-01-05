function PBM_plot(T, Y, params)

    %Recast
    mesh = params.mesh;
    mmesh = params.mmesh;

    %Determine if Y is for a single timestep, or more
    [Ny, Nt] = size(Y);

    %Interpolate result to vertices
    if params.Nr > 1 %z and r direction

        %Show final value
        for ix = 1:params.Nms
            xind = 3;
            inds = xind:params.N_rs:((params.N_cells-1) * params.N_rs + xind);
            Yf = Y(end,inds);
            Yf = reshape(Yf, Nz, Nr);
        
            if Nr == 1 %1D
                rsc = 1;
            else %2D
                rsc = mesh.volcell_cents(:,1);
                rsc = reshape(rsc, Nz, Nr);
            end
            zsc = mesh.volcell_cents(:,2);
            zsc = reshape(zsc, params.Nz, params.Nr);
            rs = unique(round(rsc(:),4));
            zs = unique(round(zsc(:),4));
            [rsc,zsc] = meshgrid(rs, zs);
            Yfi = interp2(rsc, zsc, Yf, rr, zz, 'spline');
    
            figure()
            contourf(rr, zz, Yfi);
            colorbar(); clim([0,1])
            axis equal;
        end

            
    else %Only z direction

        %Plot front
        figure();
        for ix = 1:params.Nms
            xind = ix;
            inds = xind:params.Nms:((params.N_cells-1) * params.Nms + xind);
            zsc = mesh.volcell_cents(:,2);
            zsc = reshape(zsc, params.Nz, params.Nr);
            subplot(1,2,1);
            if Nt > 1
                Yf = Y(end,inds);
            else
                Yf = Y(inds);
            end
            plot(zsc, Yf, 'LineWidth', 1.5); hold on;
        end

        %Plot aesthetics
        axis square; grid on; grid minor;
        xlabel('Z'); ylabel('Numeric Density (1/m^3)');

        %Plot bubble size distribution
        subplot(1,2,2);
        if Nt > 1
            Yf = Y(end,:);
        else
            Yf = Y;
        end
        zsc = mesh.volcell_cents(:,2);
        zsc = repmat(zsc, 1, params.Nms);
        xsc = repmat(mmesh.xsc, mesh.N_cells,1);
        Yf = reshape(Yf,params.Nms, mesh.N_cells); Yf = Yf';

        surf(100.*zsc, 1000.*xsc, Yf);
        ylabel('Bubble Mass (\mug)')
        xlabel('Vertical Position (cm)')
        zlabel('Numerical Density');
        set(gca, 'YScale', 'log');
        axis square; grid on; grid minor;

        %Create new plot showing profile at each layer
        

        %Animation plot
        if params.output.animate & Nt > 1
            figure('Units', 'Normalized', 'OuterPosition', [0, 0, 1, 1]);
            for it = 1:length(T)
    
                t = T(it);
                Yit = Y(it,:);
                sc = mesh.volcell_cents(:,2);
                zsc = mesh.volcell_cents(:,2);
                zsc = reshape(zsc, params.Nz, params.Nr);
                zsc = repmat(zsc, 1, params.Nms);
                xsc = repmat(mmesh.xsc, mesh.N_cells,1);
                Yit = reshape(Yit,params.Nms, mesh.N_cells); Yit = Yit';
    
                surf(100.*zsc, 1000.*xsc, Yit);
                title(sprintf('t = %0.4f s', t));
                ylabel('Bubble Mass (\mug)')
                xlabel('Vertical Position (cm)')
                zlabel('Numerical Density');
                axis square; grid on; grid minor;
                view([45, 45])

                %Save result
                if params.output.saveanimation

                    % Capture the plot as an image 
                    frame = getframe(gcf); 
                    im = frame2im(frame); 
                    [imind,cm] = rgb2ind(im,256); 
                    
                    % Write to the GIF File 
                    if it == 1 
                      imwrite(imind,cm,fname,'gif', 'Loopcount',inf, 'DelayTime',params.output.gifdelay); 
                    elseif mod(it,10) == 0
                      imwrite(imind,cm,fname,'gif','WriteMode','append', 'DelayTime',params.output.gifdelay); 
                    end 
                    
                end
    
                pause(0.01);


            end
        end

        %

        x=1;
    end








end