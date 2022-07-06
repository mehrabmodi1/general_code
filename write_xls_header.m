function [c] = write_xls_header(header_cells, data_cols, save_path)
%Syntax: function [cell_array] = write_xls_header(header_cells, data_cols, save_path)
%This function writes an excel file with a header row at the top and data
%columns beneath the header row. File name and path are specified in
%save_path.

%testing lines
% save_path = 'C:\Backup\Stuff\Janelia\paper_drafts\Mehrab_papers\PaBaEl2\fig_data\MBON_transitions_fig\testdata.xlsx';
% header = {'name', 'Age'};
% data_cols(:, 1) = [1;5.02;-98;63.9;0;-.2;56];
% data_cols(:, 2) = [2;5.02;-98;63.9;0;-.2;56];

n_rows = size(data_cols, 1);
n_cols = size(data_cols, 2);
n_header_labels = size(header_cells, 2);

if isempty(header_cells) == 0
    if n_cols ~= n_header_labels
        error('Number of labels in header not equal to number of data columns.')
    else
    end

    c = cell( (n_rows + 1), n_cols);
    c(1,:) = header_cells;
    c(2:end,:) = num2cell(data_cols);
else
    c = num2cell(data_cols);
end


if save_path(end) == 's'
    xlswrite(save_path,c);
elseif save_path(end) == 'x'
    writecell(c,save_path);
else
end
    

