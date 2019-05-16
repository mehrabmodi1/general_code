function [path_out] = manage_base_paths(path_in, out_type)
%This function takes a path (typically from a path list .xls file) and
%replaces the early part (computer specific) with the appropriate, locally
%saved base path, for raw (out_type 1) or analysed data (out_type 2).


%reading in base_paths specific to this computer (convention, raw data in row1col1 and an data in row2col1)
saved_local_paths = [{'C:\Data\Data\Raw_data'};...
                     {'C:\Data\Data\Analysed_data\Manual_ROI_results'};...
                     {'C:\Data\Data\Analysed_data\Suite2P_results'}...
                     ];

if out_type == 1
    base_path = saved_local_paths{1, 1};
elseif out_type == 2
    base_path = saved_local_paths{2, 1};
elseif out_type == 3
    base_path = saved_local_paths{3, 1};
end

%parsing path_in to remove base region
path_speci = findstr(path_in, '\20');
path_spec = [path_in(path_speci:end), '\'];
path_out = [base_path, path_spec];