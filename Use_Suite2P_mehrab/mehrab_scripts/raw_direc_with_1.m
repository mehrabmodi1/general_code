function [raw_direc] = raw_direc_with_1(raw_direc_base, raw_direc)
%This function checks if the existing raw_direc path has the needed \1 folder at the end and adds it if needed.

is1i = findstr(raw_direc, '\1');

if isempty(is1i) == 1
    %checking if \1 folder exists
    add_needed = isfolder([raw_direc_base, raw_direc, '\1']);
    if add_needed == 1
        raw_direc = [raw_direc, '\1'];
    else
    end
    
else
end



