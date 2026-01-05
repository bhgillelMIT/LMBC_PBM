function Cd = CalcCd(ub, Db, rho_L, mu_L, sigma, Eob)

    %Calculate reynolds number
    if ub > 0
        Reb = (rho_L * ub * Db)/mu_L;
    else
        Reb = 1E-6;
    end

    %Calculate drag coefficient
    Cd_sph = (24/Reb)*(1 + 0.15*Reb^(0.687));
    Cd_ell_cap = 8/3 * Eob/(Eob + 4);
    Cd = max([Cd_sph, Cd_ell_cap]);

end