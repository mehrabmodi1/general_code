function [col_n] = find_stim_mat_simple_col(col_name, column_heads)
%syntax:[col_n] = find_stim_mat_simple_col(col_name)
%This function takes a column title as an input and identifies it's
%location in the column_heads vector of column names. 
found_name = 0;
for name_n = 1:length(column_heads)
    curr_name = column_heads{name_n};
    
    if strcmp(col_name, curr_name) == 1
        col_n = name_n;
        found_name = 1;
        break
    else
    end
end
if found_name == 0
    error('Column name not found - re-check.')
else
end