clear all
close all

olf_n = 1;      %olfactometer number whose PID traces you want to plot
odor_n = 1;     %odor number whose PID traces you want to plot

traces_to_check = 2;    %1 - PID, 2 - LED
start_date = [];  %script scans each directory in each dir list after specified date. If empty, then it starts with first directory in each list.

%NOTE: PA on olf2 stopped being delivered on 20201120 due to a loose
%connection. This was noticed and fixed on 20201203.

dataset_list_paths = [...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\MBONG2_PaBaEl_simp_pairing_Berry.xls'};...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\MBONG2_PaBaEl_handover_pairing_Berry.xls'};...
                      {'C:\Data\Code\general_code_old\data_folder_lists\Janelia\MBONG2_PaBaEl_handover_pairing_Berry_15s_ipi.xls'};...
                      {'C:\Data\Code\general_code_old\data_folder_lists\Janelia\MBONG2_PaBaEl_handover_pairing_Berry_ELsecond.xls'};...
                      {'C:\Data\Code\general_code_old\data_folder_lists\Janelia\MBONG2_PaBaEl_handover_pairing_Berry_longpulse2.xls'};...
                      {'C:\Data\Code\general_code_old\data_folder_lists\Janelia\Berry_handover_MB298B_MBONG4-G1G2_GcaMP6f_starved.xls'};...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\Berry_handover_13F02_gcaMP6f_starved.xls'};...
                      {'C:\Data\Code\general_code_old\data_folder_lists\Janelia\d5HT1b_similar_od_handovers.xls'};...
                      
                      ];
            
suppress_plots = 1;
plotting_quant_no_filt = 0;     %1 - only unfiltered traces used for all analysis and plotting - traces included. 0 - filtered traces used for everything.

[del, odor_names1] = xlsread('C:\Data\Code\general_code_old\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList.xls', 1);
[del, odor_names2] = xlsread('C:\Data\Code\general_code_old\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList_olf2.xls', 1);
odor_names2{3} = 'Butyl acetate';

global color_vec;                
a = colormap('bone');
global greymap
greymap = flipud(a);
fly_n = 0;
script_name = mfilename;

y_ax_traces = 0.8;
y_ax_fit_traces = 0.6;

an_save_path = 'C:\Data\Data\Analysed_data\Analysis_results\PaBaEl_Gamma2\';
force_resave = 1;

n_sec = 2;      %width of moving integration window in s

integ_win_s = 5;        %width of pulse integration window for response quantification

for list_n = 1:size(dataset_list_paths, 1)
    
    curr_dir_list_path = dataset_list_paths{list_n, 1};
    [del, dir_list] = xlsread(curr_dir_list_path, 1);        %list of Suite2P results directories
    n_dirs = size(dir_list, 1);
    dataset_list_name = findstr(curr_dir_list_path, 'list_');
    dataset_list_name = curr_dir_list_path((dataset_list_name + 5):(end - 4));
    
    if isempty(findstr(curr_dir_list_path, 'handover')) == 0
        if isempty(findstr(curr_dir_list_path, 'ELfirst')) == 0
            set_list_type = 2;  %handover stimulus with EL on pulse1 case
        elseif isempty(findstr(curr_dir_list_path, 'ELsecond')) == 0
            set_list_type = 3;  %handover stimulus with EL on pulse1 case
        elseif isempty(findstr(curr_dir_list_path, 'ELfirst')) == 1 && isempty(findstr(curr_dir_list_path, 'ELsecond')) == 1
            set_list_type = 1;  %handover stimulus with no EL pulses
        else
        end
    else 
        set_list_type = 0;  %simple stimulus case
    end
    
    
    saved_traces_all = [];
    saved_PID_traces_all = [];
    %loop to go through all experiment datasets listed in list file
    for dir_n = 1:n_dirs
        fly_n = fly_n + 1;
             
        curr_dir = [dir_list{dir_n, 1}, '\'];
        curr_dir = manage_base_paths(curr_dir, 2)
        
        curr_datei = findstr(curr_dir, '2020');
        curr_date = str2num(curr_dir(curr_datei:(curr_datei + 7) ));
        if isempty(start_date) == 0
            if curr_date < start_date
                continue        %skipping current directory if it's earlier than specified start_date
            else
            end
        else
        end
           
        
        
        tif_times = load([curr_dir, 'tif_time_stamps.mat']);           %reading in time stamps for each tif file recorded by raw_data_extracter
        tif_times = tif_times.time_stamps;
        
        cd(curr_dir);
        tif_name = dir('*.tif');
        if isempty(tif_name) == 1
            continue
        else
        end
        stack_obj = ScanImageTiffReader([curr_dir, tif_name(1).name]);
        [frame_time, zoom, n_chans, PMT_offsets] = SI_tif_info(stack_obj);
        
        
        [stim_mat, stim_mat_simple, column_heads, color_vec, good_tr_list, stim_mat_orig, PID_traces, PID_traces_orig] = load_params_trains_modular(curr_dir, tif_times, frame_time);    %reading in trial stimulus parameters after matching time stamps to F traces
        PID_traces = PID_traces(:, :, traces_to_check);
        PID_traces_orig = PID_traces_orig(:, :, traces_to_check);
                
        %identifying stim_mat_simple col numbers
        led_on_col_n = find_stim_mat_simple_col('led_on', column_heads);            %identifying relevant column number in stim_mat_simple
        od_olf1_col_n = find_stim_mat_simple_col('odor_n', column_heads);           %identifying relevant column number in stim_mat_simple
        od_olf2_col_n = find_stim_mat_simple_col('odour_olf2', column_heads);       %identifying relevant column number in stim_mat_simple
        dur_olf1_col_n = find_stim_mat_simple_col('duration', column_heads);        %identifying relevant column number in stim_mat_simple
        dur_olf2_col_n = find_stim_mat_simple_col('duration_olf2', column_heads);   %identifying relevant column number in stim_mat_simple
        od_col_ns = [od_olf1_col_n, od_olf2_col_n];
        dur_col_ns = [dur_olf1_col_n, dur_olf2_col_n];
        od_durs = unique(stim_mat_simple(:, dur_col_ns(2)));
        od_durs(isnan(od_durs)) = [];
        odn_list_olf2 = unique(stim_mat_simple(:, od_col_ns(2)));
        odn_list_olf2(isnan(odn_list_olf2)) = [];
        
        
        if traces_to_check == 1
            curr_trs = find(stim_mat_simple(:, od_col_ns(olf_n)) == odor_n);
            curr_PID_traces = PID_traces(:, curr_trs);
        elseif traces_to_check == 2
            stim_mat_simple_orig = generate_stim_mat_simple(stim_mat_orig);
            curr_trs = find(stim_mat_simple_orig(:, 18) == 1);   %since we're looking at stim LED traces
            curr_PID_traces = PID_traces_orig(:, curr_trs);           %actually LED traces
        else
        end
        
        if isempty(curr_trs) == 1
            continue
        else
        end
        
        stim_frs = compute_stim_frs_modular(stim_mat, curr_trs(1), frame_time);
        stim_frs = stim_frs{olf_n};    
        
        
        plot(curr_PID_traces, 'Color', [0.6, 0.6, 0.6]);
        set_xlabels_time(1, frame_time, 10);
        ylabel('PID signal (V)')
        add_stim_bar(1, stim_frs, [0, 0, 0]);
        
        keyboard
        close figure 1
    end
end