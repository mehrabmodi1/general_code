function [matched_celln] = lookup_cell_vec(cell_contents, cell_vector)
%syntax:[cell_n] = lookup_cell_vec(cell_contents, cell_vector)
%This function takes a cell vector as an input and identifies the cell
%number that contains cell_contents. 
matched_celln = [];
for cell_n = 1:length(cell_vector)
    curr_contents = cell_vector{cell_n};
    
    if isa(cell_contents, 'string') == 1
        if strcmp(col_name, curr_contents) == 1
            matched_celln = [matched_celln; cell_n];
        else
        end
        
        if col_name == curr_contents
            matched_celln = [matched_celln; cell_n];
        else
        end
        
    else
    end
end