function fv_ih_max_funcs = Calc_fvhi_max_func(params)
    %Use sigmas at each spatial cell, can then be precise, rather than
    %using interpolation. Simplifies it two parameters (d, e_ih);

    warning('off','all')
    options = optimset('Display','off');

    %Update user
    fprintf('-- Solving Breakage Inputs\n')

    %Debug timing
    % if params.debug
    %     fprintf('')
    % end

    %Settings
    N_pts = 1000;

    %Calculate range of lambdas
    u_spfs = linspace(params.u_spfs(1), params.u_spfs(2), 20);
    turb_eps = u_spfs .* params.g;
    nus = repmat(params.nus(:), 1, length(turb_eps));
    turb_eps = repmat(turb_eps, params.Nz, 1);
    lambda_komogorov = ((nus.^3)./turb_eps).^0.25; %m
    lambda_min = 31.4 * min(lambda_komogorov(:));
    lambdas = linspace(lambda_min, 2.*params.rms(end));
    lambda_max =  2.*params.rms(end);

    %Bound energies
    fvs = linspace(0,0.5);
    fv_min = params.break.dfv_cap;
    fv_max = 0.5;
    d_min = 2 .* params.rms(1);
    d_max = 2 .* params.rms(end);
    sigma_max = max(params.sigmas);
    sigma_min = min(params.sigmas);
    e_is_max = (pi .* lambda_max.^3 .* sigma_max)/(6 .* fv_min.^(1/3) .* d_min);
    e_is_min = (pi .* lambda_min.^3 .* sigma_min)/(6 .* fv_max.^(1/3) .* d_max);
    e_is = logspace(log10(e_is_min), log10(e_is_max), N_pts);

    %Create matrix for 2d interpolation
    ds = 2 .* params.rms;
    [ds_mat, e_is_mat] = meshgrid(ds, e_is);

    %Solve
    fv_ih_max_funcs = cell(1, params.Nz);
    for iz = 1:params.Nz
        fvih_maxes = zeros(size(ds_mat));
        for ie = 1:length(e_is)
            for id = 1:length(ds)
                e_i = e_is_mat(ie, id);
                d = ds_mat(ie, id);
                fv_ih_max_func = @(fv) e_i - (fv.^(2/3) + (1 - fv).^(2/3) - 1) .* pi .* d.^2 .* params.sigmas(iz);
                try
                    fv_ih_max = fzero(fv_ih_max_func, [0, 0.5], options);
                catch
                    fv_ih_max = 0;
                end
                fvih_maxes(ie, id) = fv_ih_max;
            end

        end

        %Store interpolation function for each z-level
        fv_ih_max_func = @(d, ie) interp2(ds_mat, e_is_mat, fvih_maxes, d, ie, 'cubic');
        fv_ih_max_funcs{iz} = fv_ih_max_func;

    end

    %Turn warnings back on
    warning('on','all')

    %Save result
    if params.break.fvih_saveoutput

        %Define output name
        Nz = params.Nz;
        TLtop = params.T_Lz(params.zms(1));
        TLbot = params.T_Lz(params.zms(end));

        outname = sprintf('fv_ih_max_funcs_Nz-%d_TLtop-%d_TLbot-%d.mat', Nz, round(TLtop), round(TLbot));

        outpath = [params.folders.break, outname];
        save(outpath, 'fv_ih_max_funcs');
    end

end


