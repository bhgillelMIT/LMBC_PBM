function PBM_saveresults(T, Y, params)

    %Define filename
    date_str = string(datetime('now'));
    date_str = strrep(date_str, ':', '-'); 
    date_str = strrep(date_str, ' ', '_');
    filename = sprintf('PBM_output_%s.mat', date_str);
    filepath = [params.output.folder, filename];

    %Define output structure
    output.T = T;
    output.Y = Y;
    output.params = params;

    %Save output
    save(filepath, 'output');

end