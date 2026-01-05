function figout = PlotMovingPivots(params, pivots, plotconfig)

    %Plot colors
    colors.hydrogen = [0, 176/255, 240/255];
    colors.methane = [234/255, 112/255, 13/255];
    colors.navyblue = [0/255,0/255,128/255];
    colors.trueblue = [0/255,115/255,207/255];
    colors.elecblue = [125/255,249/255,255/255];
    colors.raspred = [179/255,68/255,108/255];
    colors.truered = [255/255,0/255,0/255];
    colors.rubyred = [155/255,17/255,30/255];
    colors.MITRed = [163/255, 31/255, 52/255];
    colors.MITDarkGrey = [138/255, 139/255, 140/255];
    colors.MITLightGrey = [194/255, 192/255, 191/255];
    colors.H2blue = [0/255, 176/255, 240/255];
    colors.ColumbiaBlue = [155/255, 221/255, 255/255];
    colors.ChiliRed = [194/255, 24/255, 7/255];
    
    %Pull values
    lw = plotconfig.lw;
    ms = plotconfig.ms;
    fs = plotconfig.fs;
    
    %Create figure
    figout = figure();

    %Pull values
    Tbs = pivots.Y(1, 1:params.N_Ts);
    Tms = pivots.Tms(1,:);
    Fms = pivots.Fms(1,:);

    %Identify indices 
    xs_in = linspace(min(pivots.Tms(1,:)), max(pivots.Tms(1,:))); 

    icore = (params.trail.N_pivots + 1):1:(length(Fms) - params.lead.N_pivots + 1);
    itrail = 1:(params.trail.N_pivots +1);
    ilead = icore(end):length(Fms);
    Tms_spline = Tms(icore);
    Fms_spline = Fms(icore);

    s2 = csape(Tms_spline, [0, Fms_spline, 0], 'clamped');
    l_trail = @(x) interp1(Tms(itrail), Fms(itrail), x, 'linear', 'extrap');
    l_lead = @(x) interp1(Tms(icore(end):length(Fms)), Fms(icore(end):length(Fms)), x, 'linear', 'extrap');
    Fbs = ppval(s2, Tms_spline);

    %Plot core interpolation
    fnplt(s2,'k-', [min(Tms_spline), max(Tms_spline)]); hold on;

    %Plot trail interpolation
    if params.trail.N_pivots > 0
        plot(Tms(itrail), l_trail(Tms(itrail)), 'k-', 'LineWidth', lw);
        stem(Tbs(itrail), l_trail(Tbs(itrail)), 'b.', 'MarkerSize', ms+6, 'LineWidth', lw, 'Color',  colors.trueblue);
    end

    %Plot lead interpolation 
    if params.lead.N_pivots > 1
        plot(Tms(ilead), l_lead(Tms(ilead)), 'k-', 'LineWidth', lw)
    end

    %Plot linear interpolation for all points
    plot(Tms, Fms, 'k-', 'Marker', '.', 'MarkerSize', ms+12, 'Color', colors.truered);
    stem(Tbs(ilead), l_lead(Tbs(ilead)), 'b.', 'MarkerSize', ms+6, 'LineWidth', lw, 'Color',  colors.trueblue);
    
    stem(Tbs(icore), ppval(s2, Tbs(icore)), 'b.', 'MarkerSize', ms+6, 'LineWidth', lw, 'Color',  colors.trueblue);

    %stem(Tbs, Fbs, 'b.', 'MarkerSize', ms+6, 'LineWidth', lw, 'Color',  colors.trueblue);

    %Create legend
    h = zeros(1,3);
    h(1) = plot(NaN, NaN, 'k-', 'LineWidth', lw); hold on;
    h(2) = plot(NaN, NaN, 'r.', 'Marker', '.', 'MarkerSize', ms+12, 'Color', colors.truered);
    h(3) = stem(NaN, NaN, 'b.', 'MarkerSize', ms+6, 'LineWidth', lw, 'Color',  colors.trueblue);
    legend(h, 'Spline/linear fit', 'Center NDs', 'Bracket NDs (interp)');

    grid on; grid minor; axis square;
    xlabel('Temperature (K)'); ylabel('Numeric Density (ND) (1/m^3*K)');
    %legend('Spline fit', 'Center NDs', 'Bracket NDs - Interpolated')
    title('Moving Pivot - Initial Distribution');
    set(gca, 'FontSize', fs, 'FontWeight', 'bold');


end