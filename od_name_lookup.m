function [od_n] = od_name_lookup(name_list, od_name)
%This function accepts a cell array of odor names and a single odor name.
%It's output is the odor number of that odor name in the odor name list.

for od_n = 1:length(name_list)
    curr_name = name_list{od_n};
    if strcmp(od_name, curr_name) == 1
        break
    else
    end
    
    
end
