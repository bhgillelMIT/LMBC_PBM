function PBM_flow_rate_study()

%% Setup

    close all

    lw = 1.5;
    fs = 18;

%% Run simulations

if false


    t_end = 10;
    uspfs = 0.01:0.01:0.1;

    base_inputs = PBM_inputs();

    for i = 1:length(uspfs)
        uspf = uspfs(i);
        inputs = base_inputs;
        inputs.reactor.u_spf_orifice = uspf;
        inputs.reactor.V_dot = uspf .* inputs.reactor.Ac;
        inputs.sim.t_end = t_end;
        PBM_v3(inputs);
    end

end

%% Evaluate outputs

    %Load directory
    folder = 'Data/Solutions/Studies/uspf_2025-10-27/';
    files = dir(folder);
    files = files(3:end,:);

    %Create storage matrices
    mdists_norm = cell(1, length(files));


    %Iterate through files and store 
    for i = 1:length(files)

        %Load file
        filename = files(i).('name');
        filepath = [folder, filename];
        data = load(filepath); data = data.output;

        %Extract final state
        params = data.params;
        T = data.T;
        Y = data.Y;

        %Plot results
        results = PBM_postprocess(T, Y, params);

        %Log results
        mdists_norm{i} = results.mdists_norm(end,:);


    end
    x = 1


    %Plot the normalized distributions
    figure();

    colororder('reef');
    semilogx(1000*params.mms, mdists_norm{1}, 'LineWidth', lw); hold on;
    semilogx(1000*params.mms, mdists_norm{2}, 'LineWidth', lw);
    semilogx(1000*params.mms, mdists_norm{3}, 'LineWidth', lw);
    grid on; grid minor; axis square;
    xlabel('Mass (mg)')
    ylabel('Normalized BSD')
    title('Normalized BSD @ Exit')
    legend('u_{spf} = 1 cm/s', 'u_{spf} = 2 cm/s', 'u_{spf} = 3 cm/s', 'location', 'northwest');
    set(gca, 'FontSize', fs)

aend