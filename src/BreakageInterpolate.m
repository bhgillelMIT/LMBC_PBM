function [beta_eddy, beta_ratio, b_eddy] = BreakageInterpolate(d, turb_eps, iz, params)

    %Determine if isothermal
    if params.isothermal
        iz_ind = 1;
    else
        iz_ind = iz;
    end



    %Pull range from 
    d_range = params.break.funcs_d_range;
    eps_range = params.break.funcs_eps_range;

    %Identify interpolation function to use
    d_valid = d >= d_range(:,1) & d <= d_range(:,2);
    eps_valid = turb_eps >= eps_range(:,1) & turb_eps <= eps_range(:,2);
    inds_valid = d_valid & eps_valid;
    inds = find(inds_valid);

    %Check validity
    N_valid = length(inds);
    if N_valid < 1
        warning('Inputs outside of data range. Limiting value that exceeds.')
        if ~any(d_valid)
            if d < min(d_range(:));
                d = min(d_range(:));
            else
                d = max(d_range(:));
            end
        elseif ~any(eps_valid)
            if eps < min(eps_range(:))
                turb_eps = min(eps_range(:));
            else
                turb_eps = max(eps_range(:));
            end
        else
            error('Unknow error.');
        end
            
    elseif N_valid > 1
        warning('Multiple windows valid, using first one.');
        inds = inds(1);
    end

    %Calculate interpolation value
    switch params.break.model
        case 'Wang_2005' %'Wang_2005'

            %Calculate functions
            beta_eddy = params.break.funcs{iz_ind}{inds}.betas(d, turb_eps, params.fvs_norm_all);
            beta_eddy = beta_eddy(:);
            beta_ratio = params.break.funcs{iz_ind}{inds}.beta_ratio(d, turb_eps);
            b_eddy = 0;


            if ~any(inds_valid)
                warning('Inputs outside of data range. Setting breakage rate to zero and beta to 1 (uniform distribution).');
                beta_eddy = ones(size(params.fvs_norm_all));
                beta_ratio = 0;
                return;
            end 
       
        case 'Luo_Svendson_1996'

            %Calculate functions
            beta_eddy = params.break.funcs{iz_ind}{inds}.betas(d, turb_eps, params.fvs_norm_all);
            beta_eddy = beta_eddy(:);
            b_eddy = params.break.funcs{iz_ind}{inds}.b_eddy(d, turb_eps);
            beta_ratio = 0;

            if ~any(inds_valid)
                warning('Inputs outside of data range. Setting breakage rate to zero and beta to 1 (uniform distribution).');
                beta_eddy = ones(size(params.fvs_norm_all));
                beta_ratio = 0;
                return;
            end
        otherwise
            error('Breakage model not recognized.');
    end


    

    %Correct 0 value at fv = 0.5
    ind_mid = floor(length(beta_eddy)/2) + 1;
    beta_eddy(ind_mid) = beta_eddy(ind_mid-1);

    %Correct nonzero values at edges
    if strcmp(params.break.model,'Wang_2005')
        beta_eddy = beta_eddy - beta_eddy(1);
    end
    
    %Handle case of no breakage
    if all(beta_eddy < 1E-4)
        beta_eddy = ones(size(beta_eddy));
    end

end