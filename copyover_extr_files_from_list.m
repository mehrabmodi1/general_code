%This script deletes the extracted raw data mat file alone from extracted
%dataset folders at sink (laptop/server) as per the lists in list_of_lists. It then copies
%over only the extracted raw data mat file from source(server/an_machine).

transfer_type = 1;      %valid options are 1 - source:an_machine, sink:server; or 2 - source:server, sink:laptop

if transfer_type == 1
    source_path_base = 'E:\Data\Analysed_data\Manual_ROIs\';    %an_machine path
    sink_path_base = '\\dm11\turnerlab\Mehrab\Analysed_data\';  %server path
    list_path_base = 'E:\Data\Raw_Data_Current\dataset_lists\';
elseif transfer_type == 2
    source_path_base = '\\dm11\turnerlab\Mehrab\Analysed_data\';  %server path
    sink_path_base =  'C:\Data\Data\Analysed_data\Manual_ROI_results\'; %laptop path
    list_path_base = 'C:\Data\Code\general_code_old\data_folder_lists\Janelia\';
else
end
 


%dataset list files listed here according to laptop convention.
list_of_lists =  [{'dataset_list_Alpha1_60strace_71C03LxA_MB043CGal4_Chrison_noLED_control.xls'};...
    {'dataset_list_Alpha1_60strace_71C03LxA_MB043CGal4_noChrisoncontrol.xls'};...
    {'dataset_list_Alpha1_60strace_71C03LxA_MB043CGal4_with_Chrimson_LED.xls'};...
    {'dataset_list_MBONG2_PaBaEl_handover_starved_halfAra.xls'};...
    {'dataset_list_MBONG2_PaBaEl_handover_starved_halfAra_prehabituated.xls'};...
    {'dataset_list_MBONG2_PaBaEl_handover_starved_halfAra_prehabituated_strongUS.xls'};...
    {'dataset_list_MBONG2_PaBaEl_handover_starved_halfAra_prehabituated_strongUS_EL_handover.xls'};...
    {'dataset_list_MBONG2_PaBaEl_handover_starved_halfAra_prehabituated_strongUS_EL_second.xls'}];


for list_n = 1:size(list_of_lists, 1)
    curr_list_name = list_of_lists{list_n, 1};
    if transfer_type == 1
        curr_list_name = curr_list_name(14:end);
    else
    end
    
    curr_list_path = [list_path_base, curr_list_name];
    [del, dir_list] = xlsread(curr_list_path, 1);        %list of results directories
    n_dirs = size(dir_list, 1);
    
    for dir_n = 1:n_dirs
        curr_dir = [dir_list{dir_n, 1}, '\'];
        curr_dir_orig = curr_dir;
        curr_dir_slashi = findstr(curr_dir_orig, '\');
        curr_dir = curr_dir(curr_dir_slashi(end-2):curr_dir_slashi(end));
        
        sink_path = [sink_path_base, curr_dir, 'extracted_raw_data_mat.mat'];
        delete(sink_path);
        
        source_path = [source_path_base, curr_dir, 'extracted_raw_data_mat.mat'];
        
        %copying over new extracted_raw_data_mat.mat file for curr_dir from file server
        status = copyfile(source_path, sink_path);
        disp([num2str(status), '- copy status for ', curr_dir]);
       
    end
    
    
end