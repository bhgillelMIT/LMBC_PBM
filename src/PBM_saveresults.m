function PBM_saveresults(T, Y, params)

    %Define filename
    date_str = string(datetime('now'));
    date_str = strrep(date_str, ':', '-'); 
    date_str = strrep(date_str, ' ', '_');
    temp_react_str = sprintf('T=%d_X=%d', double(params.heat.active),double(params.react.active));
   
    filename = sprintf('PBM_%s_%s.mat', temp_react_str, date_str);
    filepath = [params.output.folder, filename];

    %Define output structure
    output.T = T;
    output.Y = Y;
    output.params = params;

    %Save output
    save(filepath, 'output');

end