%User params
curr_dir = 'C:\Data\Data\Analysed_data\Suite2P_results\20180628\fly3_c739AB_x_opGC6f_Yoshi_PaBaEl\1\';



[del, odor_names] = xlsread('C:\Data\Code\general_code\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList.xls', 1);
global color_vec;                
color_vec = load('C:\Users\Mehrab\Google Drive\Backup\Stuff\CSHL\Glenn lab\Code\std_color_vec.txt');

%odors_n;  if empty, all odors will be plotted one by one
odor_list = [];

%dur;   if empty, all durations will be plotted one by one
odor_dur_list = 60;

[stim_mat, stim_mat_simple, column_heads, color_vec] = load_params_trains(curr_dir, []);    %reading in trial stimulus parameters after matching time stamps to F traces
       
%Reading in experimental parameters
if isempty(odor_list) == 1
    odor_list = unique(stim_mat_simple(:, 2) );
else
end

n_odors = length(odor_list);

if isempty(odor_dur_list) == 1
    odor_dur_list = unique(stim_mat_simple(:, 3) );
else
end

n_od_durs = length(odor_dur_list);
n_trains = max(stim_mat_simple(:, 11));
saved_an_results.odor_list = odor_list;
saved_an_results.odor_dur_list = odor_dur_list;

cd(curr_dir);
tif_name = dir('*.tif');
stack_obj = ScanImageTiffReader([curr_dir, tif_name(1).name]);
[frame_time, zoom, n_chans, PMT_offsets] = SI_tif_info(stack_obj);


%loading extracted raw fluorescence data matrices written by raw_dff_extractor
raw_data_mat = load([curr_dir 'extracted_raw_data_mat.mat']);
raw_data_mat = raw_data_mat.raw_data_mat;           %raw F traces extracted from ROIs
raw_data_mat_orig = raw_data_mat;
raw_data_mat = raw_data_mat(:, :, stim_mat_simple(:, 1));       %making sure only time-stamp matched trials are used for further analysis
n_cells = size(raw_data_mat, 2);

%calculating dF/F traces from raw data
filt_time = 0.5;            %in s, the time window for boxcar filter for generating filtered traces
[dff_data_mat, dff_data_mat_f] = cal_dff_traces_res(raw_data_mat, stim_mat, frame_time, filt_time, curr_dir);
del = find(dff_data_mat_f < -1);
dff_data_mat_f(del) = -1;       %forcing crazy values to sane ones

%identifying significantly responsive cells
[resp_sizes, sig_trace_mat, sig_trace_mat_old, sig_cell_mat] = cal_sig_responses_res(dff_data_mat, stim_mat, stim_mat_simple, curr_dir, frame_time);
sig_cells = find(sum(sum(sig_cell_mat, 3), 2) > 0);         %list of all cells significant for any odor for any duration


%plotting response traces
for odor_n = 1:length(odor_list)
    odor_ni = odor_list(odor_n);
    true_odor_list = unique(stim_mat_simple(:, 2));
        true_odor_n = find(true_odor_list == odor_ni);     
    for odor_dur_n = 1:length(odor_dur_list)
        curr_odor_dur = odor_dur_list(odor_dur_n);
        true_dur_list = unique(stim_mat_simple(:, 3));
        true_dur_n = find(true_dur_list == curr_odor_dur);
        
        curr_sig_cells = find(sig_cell_mat(:, odor_ni, true_dur_n) == 1); 
        for cell_n = 1:length(curr_sig_cells)
            cell_ni = curr_sig_cells(cell_n);
            [h, max_val] = plot_traces_simple(1, dff_data_mat, cell_ni, true_odor_n, true_dur_n, stim_mat_simple, stim_mat, frame_time, color_vec, 0, 1);
            keyboard
            close figure 1
            
        end
        
    end
end
    
