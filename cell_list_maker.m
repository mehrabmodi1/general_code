function [cell_list] = cell_list_maker(direc, odor_ni)

if exist([direc 'cell_list_odor_' int2str(odor_ni) '.txt']) == 2                   %pooling interesting cells for a given odor across different stim durations
    cell_list = load([direc 'cell_list_odor_' int2str(odor_ni) '.txt']);
else
    cell_list = [];
end
end_loop = 0;
%loop to let user pick interesting cells
while end_loop == 0
    figure(1)
    cell_n = input('Enter cell number to keep track of. To stop, enter 0.');
    if cell_n ~= 0
        cell_list = [cell_list; cell_n];
    else
        end_loop = 1;
    end
end
    
cell_list = unique(cell_list);      %making sure repeat runs of program don't cause duplicate entries

save([direc 'cell_list_odor_' int2str(odor_ni) '.txt'], 'cell_list', '-ASCII');


end