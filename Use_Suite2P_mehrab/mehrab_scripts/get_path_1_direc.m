%checking if directory structure was extended with a \1\ folder at the
%end and copying over the results file from ROI_prune to that folder

function [path_out] = get_path_1_direc(path_in)
curr_path = path_in;
end_reached = 0 ;

%in case the folder doesn't exist yet, there are no \1 direcs and path_in
%is returned
if isdir(curr_path) == 1

    %loop to get to last \1 direc
    while end_reached == 0
        cd(curr_path);
        if isdir([curr_path, '1\']) == 1
            curr_path = [curr_path, '1\'];
        else
            end_reached = 1;
        end
    end
    
else
end

path_out = curr_path;