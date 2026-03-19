% A function implementing the turbulent eddy coalescence rate model of Hibiki and Ishii 2002. 
% Available here: https://www.sciencedirect.com/science/article/pii/S0017931001003271?via%3Dihub

function Coalescence_HibikiIshii(y, j,k, zind, params, turb)

    %Define index to sample from for liquid properties
    i = (zind - 1) * params.Nms + 1;

    %Pull mass values at pivots
    mj = params.mms(j);
    mk = params.mms(k);

    %Determine which bracket 

    %Pull other values at pivots
    jind = j; %jind = (zind - 1) .* params.Nms + j;
    kind = k; %(zind - 1) .* params.Nms + k;
    Nj = y(jind);
    Nk = y(kind);

    %Calculate coalesced mass
    mjk = mj + mk;

    %Calculate diameters of bubbles 
    dj = params.d_mu(zind, jind); %((6 * mj * params.R .* params.Ts(i))./(pi .* params.ps(i) .* params.Ms(i))).^(1/3);
    dk = params.d_mu(zind, kind); %((6 * mk * params.R .* params.Ts(i))./(pi .* params.ps(i) .* params.Ms(i))).^(1/3);

    %Calculate collision frequency based on Hibiki and Ishii 2002
    Lo = sqrt(params.sigmas(zind)./(params.g .* params.rhos_l(zind)));
    Re = ((turb.eps(zind).^(1/3) .* Lo.^(1/3)).*Lo)/params.nus(zind);
    gamma_c = (1.82E-8) .* params.alpha_g(zind) .* Re.^3;
    omega_c = params.alpha_g(zind) .* turb.eps(zind).^(1/3)
    
    %Calculate coalescence probability based on Hibiki and Ishii 2002
    Kc = 1.29; %Model parameter - Hibiki and Ishii 2002
    P_c = exp(-(Kc .* params.rho_l(zind).*)./());



    %Calculate coalescence rate
    c_jk = P_c .* omega_c;

end