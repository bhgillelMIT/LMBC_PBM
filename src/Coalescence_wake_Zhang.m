%Model derived from "A theoretical model for coalescence efficiency in
%collisions..." by Xi-Bao Zhang (2025)

function cw_jk = Coalescence_wake_Zhang(dj, dk, turb, zind, i, params)

    %Pull
    dlead = max([dj, dk]); %m - diameter of the larger bubble
    dfollow = min([dj, dk]);


    %Calculate contact time
    theta = 50/pi;
    Ri = dlead/(2*(1-cos(theta))^2 * (2+cos(theta))).^(1/3);
    beta = 0.25;
    eta = 0.9;
    Rc = dfollow/2; %m - radius of following bubble
    rc = Rc * sin(theta); %m - radius of contact 
    u_ro = 
    

    t_contact = (1./(2 * beta * Rc*u_ro)) .* (eta*m_eq*u_ro^2./(4*params.sigmas(iz)*pi) + rc^2);


    %Calculate drainage time
    t_drain = sqrt(rho_l * rc^2./(4*sigma*(1/Rc + 1/rc))) * log(2*h0/hc);


    %Calculate coalescence efficiency
    Pw_jk = exp(-t_drain./t_contact);

    %Calculate collision frequency
    K_w1 = 1.3;
    dc_w = 4 .* sqrt(params.sigmas(zind)./(params.g .* params.rhos_l(i))); %Critica diameter for wake entrainment
    if dlead > dc_w/2

        %Calculate slip velocity - seems to be for a cap bubble
        u_bar_slip = 0.71 .* sqrt(params.g.*dk);
    
        %Determine if bubble is above critical size for wake
        %entrainment
        if dlead >= dc_w/2
            theta_c = ((dlead - dc_w/2).^6)./((dlead - dc_w/2).^6 + (dc_w/2).^6);
        else
            theta_c = 0;
        end
    
        %Calculate coalescence rate
        cw_jk = K_w1 .* theta_c .* dlead.^2 .* u_bar_slip .* Pw_jk;
    else
        cw_jk = 0;
    end




end