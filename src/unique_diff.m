function [vals, inds] = unique_diff(A, thresh)

    

    %Iterate through and take only values separated by more than thresh
    vals = A(1);
    inds = 1;
    for ia = 2:length(A)
        val = A(ia);
        diff = abs(val - vals(end));
        if diff > thresh
            vals(end+1) = val;
            inds(end+1) = ia;
        end

    end

        


    


end