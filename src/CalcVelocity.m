function ub = CalcVelocity(Db, rho_L, mu_L, sigma, Eob, fsolve_opts)

    

    %Physical constants
    g = 9.81;

    %Initial velocity guess
    ub_guess = 0.3;

    ub_func = @(ub) sqrt((4*Db*g)./(3 * CalcCd(ub, Db, rho_L, mu_L, sigma, Eob))) - ub;
    ub = fsolve(ub_func, ub_guess, fsolve_opts);
    
    


end