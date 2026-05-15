function params = LoadBreakFile(params)

    %global params

    if params.break.loadfile
        try

            %Load file
            load(params.break.filename);

            %Pull ranges for interpolation functions
            params.break.funcs_d_range = break_out.d_range;
            params.break.funcs_eps_range = break_out.eps_range;

            %Check if the number of z-layers matches
            Nz_file = size(break_out.funcs, 1); %Pull the number of z-layers in the breakage file
            if Nz_file ~= params.Nz & params.break.stretch_model
                for iz = 1:params.Nz

                    %Pull breakage functions for this z-layer
                    if Nz_file == 1
                        break_out.funcs{iz} = break_out.funcs{1}; %If only one z-layer in file, use for all layers
                    else
                        error('Interpolation between breakage layers not yet implemented; Use single layer file instead.'); %break_out.funcs{iz} = break_out.funcs{iz}; %Use corresponding z-layer if number of layers matches
                    end
                    
                end
                params.break.funcs = break_out.funcs;
            elseif Nz_file ~= params.Nz & ~params.break.stretch_model
                error("Breakage file does not have the same number of z-layers as the current simulation. Check file or set break.loadfile to false.");
            else
                params.break.funcs = break_out.funcs;
                
            end




            


            %params.break.beta_ratio = break_out.beta_ratio;
            %params.break.beta = break_out.beta;
        catch
            params = CalculateBSDs(params);
        end
    else
        params = CalculateBSDs(params);
    end


end