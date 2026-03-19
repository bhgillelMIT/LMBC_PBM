function SourceTest(params)

    %Settings
    debug = true;
    lw = 2;

    %Create debug figure
    if debug
        figure();
    end

    %Make params global
    global params


    %Update input values to simplify analysis
    Nz_actual = params.Nz;
    src_its_actual = 0;
    params.Nz = 1;
    params.src.its = 1;
    params.turb.k = 0.1;
    yin = params.N_dot_os./trapz(params.N_dot_os);

    %yin = ones(1, params.Nms); %Assume a distribution matching the input distribution

    

    %Define parameters for tests
    epss = 0.1:0.1:1;

    %Storage matricies
    cadds = zeros(length(epss), params.Nms);
    badds = zeros(length(epss), params.Nms);
    hs = zeros(length(epss), params.Nms);


    %Iterate through each value of eps and calculate the source term for each size/spatial cell
    for i = 1:length(epss)

        %Update parameters
        params = CalcLocalProperties(yin, params);
        params.turb.eps = epss(i);
        
        %Calculate coalescence and breakage matrices
        [cadd, csub, cmats] = Coalescence(yin, params); %Coalescence(Fs, params );  
        [badd, bsub, bmats] = Breakage(yin, params);
        h = cadd(:) + badd(:) - csub(:) - bsub(:);

        %Log results
        cadds(i,:) = cadd';
        badds(i,:) = badd';
        hs(i,:) = h';

        %Plot the results for this eps
        if debug

            %Plot individual components
            cla
            plot(params.mms, cadd, 'r-', 'LineWidth', lw); hold on;
            plot(params.mms, csub, 'r--', 'LineWidth', lw);
            plot(params.mms, badd, 'b-', 'LineWidth', lw);
            plot(params.mms, bsub, 'b--', 'LineWidth', lw);
            plot(params.mms, h, 'k-', 'LineWidth', lw);
            xlabel('Representative Mass'); ylabel('Source Term'); title(['\epsilon = ', num2str(epss(i))]);
            grid on; grid minor; axis square;
            legend('cadd', 'csub', 'badd', 'bsub', 'h_net');
        end

        

    end

    %Restore values for params

end