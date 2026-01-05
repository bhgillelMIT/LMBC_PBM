function d_crit = CalcCritDiameter(iz, params)

    d_crit_func = @(d) WeFunc(d, iz, params) - params.coalesce.We_crit;
    
    d_crit = fzero(d_crit_func, [0, 0.3]); %(params.We_crit * params.sigmas(iz))./(params.rhos(iz) .* u.^2); %Critical diameter for 


end


function We = WeFunc(d, iz, params)
    Eob = (params.rhos(iz) .* d.^2 .* params.g)/params.sigmas(iz);
    u = CalcVelocity(d, params.rhos(iz), params.mus(iz), params.sigmas(iz), Eob, params.fsolve_opts);
    We = (params.rhos(iz) .* u.^2 .* d)/params.sigmas(iz);
    
end