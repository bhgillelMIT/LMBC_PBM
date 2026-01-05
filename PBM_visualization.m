function PBM_visualization()

    %Define reactor geometry
    r_wall_o        = 0.027;
    r_wall_i        = 0.019;
    r_base          = 0.05;
    h_base          = 1 * 0.0254;
    base_coords     = [r_base, 0;
                       r_base, h_base;
                       -r_base, h_base;
                       -r_base, 0];
        
    top_coords      = base_coords;
    top_coords(:,2) = top_coords(:,2) + 11 * 0.0254;
    
    left_coords     = [-r_wall_o, 01 * 0.0254;
                       -r_wall_o, 11 * 0.0254;
                       -r_wall_i, 11 * 0.0254;
                       -r_wall_i, 01 * 0.0254];
                   
    right_coords    = left_coords;
    right_coords(:,1) = -right_coords(:,1);
    
    %Define Tin geometry
    tin_coords      = [r_wall_i, 0;
                       r_wall_i, 10*h_base;
                       -r_wall_i, 10*h_base;
                       -r_wall_i, 0];
    
    %Define air geometry
    air_coords      = [r_wall_i, 0;
                       r_wall_i, 12*h_base;
                       -r_wall_i, 12*h_base;
                       -r_wall_i, 0];
    
    %Initialize plot and display 
    figure('units','normalized','outerposition',[0 0 0.66 1]); pause(0.1);
    fill(base_coords(:,1), base_coords(:,2), 'k', 'EdgeColor', 'none'); hold on;
    fill(top_coords(:,1), top_coords(:,2), 'k', 'EdgeColor', 'none');
    fill(left_coords(:,1), left_coords(:,2), 'k', 'EdgeColor', 'none');
    fill(right_coords(:,1), right_coords(:,2), 'k', 'EdgeColor', 'none');
    fill(air_coords(:,1), air_coords(:,2), [1, 1, 1], 'EdgeColor', 'none');
    fill(tin_coords(:,1), tin_coords(:,2), [0.5, 0.5, 0.5], 'EdgeColor', 'none');
    set(gcf,'color','w');
    
    
    %Plot strata geometry
    dz = mean(diff(zs(1,:)));
    
    zbounds = 0:dz:(nzs * dz);
    for i = 1:nzs
        zbu = zbounds(i+1);
        zbl = zbounds(i);
        rectangle('Position', [-r_wall_i, zbl, 2*r_wall_i, (zbu - zbl)]);
        
        
    end
    
    %Create bar graph outlines
    bg_start = 4*r_wall_i;
    bg_indent = 0.02;
    bg_width = 4*r_wall_i;
    bg_barwidth = bg_width*(1-2*bg_indent)/nxs;
    for i = 1:nzs
        zbu = zbounds(i+1);
        zbl = zbounds(i);
        rectangle('Position', [bg_start, zbl, bg_width, (zbu - zbl)]);
        
        
    end
    
    %Add annotations
    
    
    %Create storage vector cell array for rectangles
    rects = {}
    
    
    %Iterate through times and plot x distributions in the bar graph
    %outlines
    
    it = 1;
    for tind = tinds

        %Delete all rectangles
        for r = 1:length(rects)
            rect = rects{r};
            delete(rect);
        end
        rects = {};
        
        
        %Iterate through strata
        for i = 1:nzs
            
            zbu = zbounds(i+1);
            zbl = zbounds(i);
            
            %Pull x distribution for the current stratum and time
            xdist = Fs_output(:, i, tind);
            xdist(end-1:end) = 0.25 * xdist(end-1:end);
            xmax = max(xdist);
            if xmax > 0;
                xdist = 0.9 * dz .* xdist./xmax;
            end
            
            %Define plot limits
            hmin = zbl;
            hmax = hmin + dz*0.9;
            
            %Iterate through x dist and plot in 
            xi = bg_start + bg_indent*bg_width;
            for xind = 1:nxs
                
                rects{end + 1} = rectangle('Position', [xi, zbl, 0.9*bg_barwidth, xdist(xind)], 'FaceColor', 'blue', 'EdgeColor', 'none');
                xi = xi + bg_barwidth;
            end

        end
        
        %Annotations
        ttext = text(0.1, 0.3, ['t = ', num2str(ts_output(tind)), ' s'], 'FontSize', 24, 'FontWeight', 'bold'); 
        
        %Draw
        axis off 
        axis equal
        drawnow 
        
        % Capture the plot as an image 
        frame = getframe(gca); 
        im = frame2im(frame); 
        [imind,cm] = rgb2ind(im,256); 
        % Write to the GIF File 
        if ts_output(tind) <= t_display_max
            if it == 1 
              imwrite(imind,cm,fname_2D,'gif', 'Loopcount',inf, 'DelayTime',gifdelay); 
            else 
              imwrite(imind,cm,fname_2D,'gif','WriteMode','append', 'DelayTime',gifdelay); 
            end 
        end
        

        pause(t_delay)
        it = it + 1;
        delete(ttext)
        
    end 
    
    %Update axis appearance
    
end