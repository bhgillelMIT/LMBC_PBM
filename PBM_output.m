function status = PBM_output(t, y, flag)

    %Make parameters global
    global params

    %Allocate storage vectors
    %persistent its 
    

    % if strcmp(flag, 'init')
    %     params.its = 1;
    % else
    params.its = params.its + 1;
    % end

        


    if isempty(flag)
        flag = 'None';
    elseif strcmp(flag, 'done')
        x = 1;
        %paramsits = 0; %Reset storage vector
    end


    %Calculate mass and record
    if ~strcmp(flag, 'done')
        params.m_total(params.its+1) = sum(params.mms_rep .* y(:,end) .* params.Vcells_rep);
    end
    


    %Print update
    if params.sol.solve_details
        detail_str = 'true';
    else
        detail_str = 'false';
    end
    fprintf('\nPBM (Iteration %d; flag = %s; detail = %s) \n\n', params.its, flag, detail_str);

    %Update params
    if mod(params.its, params.sol.solve_detail_its) > 0 & ~any(isnan(params.break.badd))
        params.sol.solve_details = false; %Switch off
    else
        params.sol.solve_details = true; %Switch on
    end


    %Define status
    status = 0;
end