function params = LoadBreakFile(params)

    global params

    if params.break.loadfile
        try
            load(params.break.filename);
            params.break.funcs = break_out.funcs;
            params.break.funcs_d_range = break_out.d_range;
            params.break.funcs_eps_range = break_out.eps_range;


            %params.break.beta_ratio = break_out.beta_ratio;
            %params.break.beta = break_out.beta;
        catch
            params = CalculateBSDs(params);
        end
    else
        params = CalculateBSDs(params);
    end


end