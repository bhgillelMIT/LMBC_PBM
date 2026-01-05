function SolverComparison()

    lw = 2;

    %Load data
    data1 = load('Data/Solutions/Solver_Comp_2025-12-05/PBM_output_05-Dec-2025_13-20-46.mat');
    data2 = load('Data/Solutions/Solver_Comp_2025-12-05/PBM_output_05-Dec-2025_13-44-54.mat');

    %Isolate final curves
    mms = data1.output.params.mms;
    Nm = data1.output.params.Nms;
    Y1 = data1.output.Y;
    Y2 = data2.output.Y;
    Y1_end = data1.output.Y(end,181:200)
    Y2_end = data2.output.Y(end,181:200)


    %zinds = params.zinds;


    results1 = PBM_postprocess(data1.output.T, data1.output.Y, data1.output.params)
    results2 = PBM_postprocess(data2.output.T, data2.output.Y, data2.output.params)

    %Plot
    figure();

    
    semilogx(mms.*1000, Y1_end, 'r-', 'LineWidth', lw); hold on;
    semilogx(mms.*1000, Y2_end, 'b-', 'LineWidth', lw); 
    xlabel('Bubble Mass (g)'); ylabel('Numeric Density');
    grid on; grid minor; axis square;
    title('Solver Comparison');
    legend('Direct', 'Segregated')

    %Surface plot - unnormalized
    figure()
    subplot(1,2,1);
    params = data1.output.params;
    mesh = params.mesh;
    zsc = mesh.volcell_cents(:,2);
    zsc = repmat(zsc, 1, params.Nms);
    xsc = repmat(params.mmesh.xsc, mesh.N_cells,1);
    Yf = Y1(end,:); Yf = reshape(Yf, params.Nms, params.Nz); Yf = Yf';
    %Yf = reshape(Yf,params.Nms, mesh.N_cells); Yf = Yf';
    surf(100.*zsc, 100.*xsc, Yf);
    colormap jet
    ylabel('Bubble Mass (\mug)')
    xlabel('Vertical Position (cm)')
    zlabel('Normalized BSD');
    set(gca, 'YScale', 'log');




end