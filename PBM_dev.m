function PBM_dev()

% This file is just used to write code while the main files are used for
% simulations. The code can then be copied to the file when simulations
% aren't running.

    %
    uj = params.uz_mu(iz, jind)
    uk = params.uz_mu(iz, kind)

    Pr_jk = 0.5; %Efficiency of coalescence for j and k due to different rise velocities
    cr_jk = (pi/4) .* ((dj + dk).^2) .* abs(uj - uk)*Pr_jk
    


end