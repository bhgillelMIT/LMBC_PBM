% Based off the paper "Bubble Coalescence and Break-Up in Air-Sparged
% Bubble Columns" by Prince and Blanch 1990 
% https://aiche.onlinelibrary.wiley.com/doi/pdf/10.1002/aic.690361004 

function [c_jk, c_jks] = Coalescence_PrinceBlanch(y, j,k, zind, params, turb)

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

    %Calculate turbulent eddy velocities
    u_bar_j = sqrt(2) .* (turb.eps(zind) .* dj).^(1/3); %Mean turbulent velocity 
    u_bar_k = sqrt(2) .* (turb.eps(zind) .* dk).^(1/3);

    %Calculate collision frequency based on Prince and Blanch 1990
    omega_c = (pi/4) .* (dj + dk).^2 .* sqrt(u_bar_k.^2 + u_bar_j.^2);

    %Calculate coalescence efficiency
    psi = 1;
    u_bar_jk = sqrt(u_bar_j.^2 + u_bar_k.^2);
    We_ij = (params.rhos_l(i) .* (dj) .* u_bar_jk.^2)./params.sigmas(zind);
    rho_g = 0.5;
    gamma_VM = 0.5; %coefficient of virtual mass
    Pc_jk = exp(-psi .* sqrt(0.75 .* (1+ (dj/dk).^2) .* (1 + (dj./dk).^3))./((rho_g/params.rhos_l(i) + gamma_VM) .* (1 + (dj./dk)).^3) .* sqrt(We_ij));

    %Wake entrainment
    cw_jk = 0;
       % cw_jk = 0.0073 * 

    %Add a term for fine coalescenc
    dsmaller = min([dj,dk]); dlarger = max([dj, dk]);
    cscalar = 1; %1 + 100*exp(-dlarger/params.dms(1));


    %Define outputs
    c_jk = Pc_jk .* omega_c .* cscalar + cw_jk;
    c_jks = [c_jk, 0, 0];





end