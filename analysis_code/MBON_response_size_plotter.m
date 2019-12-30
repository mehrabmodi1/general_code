clear all
close all

dataset_list_paths = [...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_fore_distr_MBONAlpha1.xls'};...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_fore_distr_MBONAlpha1_set2.xls'};...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_fore_distr_MBONAlpha1_set3_highLED.xls'};...
%                       {'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_fore_distr_MBONGamma2.xls'};...
%                       {'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_fore_distr_MBONGamma2_set1.xls'};...
%                       {'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_fore_distr_MBONGamma2_set2.xls'};...
%                       {'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_PaBaEl_MBONG2_handover.xls'};...
%                       {'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_MBONG2_PaBaEl_handover_simple_starved.xls'};...
%                       {'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_MBONG2_PaBaEl_handover_starved_set2.xls'};...
                       %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_Alpha1_60strace_71C03LxA_MB043CGal4.xls'};...
%                       {'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_fore_distr_MBONAlpha1.xls'};...
%                       {'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_fore_distr_MBONAlpha1_set2.xls'};...
%                       {'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_fore_distr_MBONAlpha1_set3_highLED.xls'};...
                        
                        {'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_MBONG2_PaBaEl_handover_starved_halfAra.xls'};...
                        %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_MBONG2_PaBaEl_handover_starved36_halfAra.xls'};...



                    ];
            
suppress_plots = 0;
[del, odor_names1] = xlsread('C:\Data\Code\general_code_old\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList.xls', 1);
[del, odor_names2] = xlsread('C:\Data\Code\general_code_old\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList_olf2.xls', 1);
odor_names2{3} = 'Butyl acetate';

global color_vec;                
a = colormap('bone');
global greymap
greymap = flipud(a);
fly_n = 0;
script_name = mfilename;

saved_data_mat = [];
n_flies = 0;
for list_n = 1:size(dataset_list_paths, 1)
    curr_dir_list_path = dataset_list_paths{list_n, 1};
    [del, dir_list] = xlsread(curr_dir_list_path, 1);        %list of Suite2P results directories
    n_dirs = size(dir_list, 1);
    dataset_list_name = findstr(curr_dir_list_path, 'list_');
    dataset_list_name = curr_dir_list_path((dataset_list_name + 5):(end - 4));
    
    saved_resps = zeros(n_dirs, 7) + nan;                    %mean response of each fly
    
    %loop to go through all experiment datasets listed in list file
    saved_resps = [];
    for dir_n = 1:n_dirs
        fly_n = fly_n + 1;
              
        saved_an_results.scriptname = mfilename('fullpath');
        curr_dir = [dir_list{dir_n, 1}, '\'];
        curr_dir = manage_base_paths(curr_dir, 2);
       
        tif_times = load([curr_dir, 'tif_time_stamps.mat']);           %reading in time stamps for each tif file recorded by raw_data_extracter
        tif_times = tif_times.time_stamps;
        try
            [stim_mat, stim_mat_simple, column_heads, color_vec] = load_params_trains_modular(curr_dir, tif_times);    %reading in trial stimulus parameters after matching time stamps to F traces
        
        catch
            [stim_mat, stim_mat_simple, column_heads, color_vec] = load_params_trains(curr_dir, tif_times);
        end
            
        %identifying stim_mat_simple col numbers
        led_on_col_n = find_stim_mat_simple_col('led_on', column_heads);            %identifying relevant column number in stim_mat_simple
        od_olf1_col_n = find_stim_mat_simple_col('odor_n', column_heads);           %identifying relevant column number in stim_mat_simple
        dur_olf1_col_n = find_stim_mat_simple_col('duration', column_heads);        %identifying relevant column number in stim_mat_simple
                
        %Reading in experimental parameters
        odor_list_olf1 = unique(stim_mat_simple(:, od_olf1_col_n) );
        
        cd(curr_dir);
%         tif_name = dir('*.tif');
%         stack_obj = ScanImageTiffReader([curr_dir, tif_name(1).name]);
%         [frame_time, zoom, n_chans, PMT_offsets] = SI_tif_info(stack_obj);
        frame_time = 0.099;

        %loading extracted raw fluorescence data matrices written by raw_dff_extractor
        raw_data_mat = load([curr_dir 'extracted_raw_data_mat.mat']);
        raw_data_mat = raw_data_mat.raw_data_mat;           %raw F traces extracted from ROIs
        raw_data_mat_orig = raw_data_mat;
        tif_n_col_n = find_stim_mat_simple_col('matched_tif_n', column_heads);
        raw_data_mat = raw_data_mat(:, :, stim_mat_simple(:, tif_n_col_n));       %making sure only time-stamp matched trials are used for further analysis
        n_cells = size(raw_data_mat, 2);
        
        %computing dF/F traces from raw data
        filt_time = 0.5;            %in s, the time window for boxcar filter for generating filtered traces
        [dff_data_mat, dff_data_mat_f] = cal_dff_traces_res(raw_data_mat, stim_mat, frame_time, filt_time, curr_dir);
        
        if size(dff_data_mat, 2) > 1
            dff_data_mat = mean(dff_data_mat, 2, 'omitnan');
            dff_data_mat_f = mean(dff_data_mat_f, 2, 'omitnan');
        else
        end
        del = find(dff_data_mat_f < -1);
        dff_data_mat_f(del) = -1;       %forcing crazy values to sane ones
        
        for tr_n = 1:size(dff_data_mat_f, 3)
            %computing response sizes
            stim_frs = compute_stim_frs_modular(stim_mat, tr_n, frame_time);
            stim_frs = stim_frs{1};      
            curr_trace = dff_data_mat_f(:, :, tr_n);
            curr_trace = squeeze(mean(curr_trace, 2, 'omitnan'));
            pk_resp = max(curr_trace(stim_frs(1):(stim_frs(2) + round(5./frame_time)) ));
            area_resp = mean(curr_trace(stim_frs(1, 1):(stim_frs(1, 2) + round(5./frame_time)) ), 'omitnan');
            
            %keeping track of total session time, imaging time
            if tr_n == 1
                tstamp1 = tif_times(tr_n).tstamp;
                tot_scan_dur = 0;
            else
            end
            curr_tstamp = tif_times(tr_n).tstamp;
            curr_time = seconds(abs(curr_tstamp - tstamp1));      %time at tr beginning
            
            tot_od_dur = sum(sum(stim_mat(tr_n).pulse_train));
            scan_dur = stim_mat(tr_n).stimLatency + tot_od_dur + stim_mat(tr_n).post_od_scan_dur;
            
            %logging data
            saved_data_mat = [saved_data_mat; pk_resp, area_resp, curr_time, tot_scan_dur];
            
            tot_scan_dur = tot_scan_dur + scan_dur;
            
           
        end
        n_flies = n_flies + 1;
    end
    
    
    

end

extra_trsi = find(saved_data_mat(:, 3) > 2000);
saved_data_mat(extra_trsi, :) = [];

%computing binned mean values and SEs
bin_vec = 0:(2000/20):2000;
for bin_n = 1:(length(bin_vec) - 1)
    curr_bini = find(saved_data_mat(:, 4) > bin_vec(bin_n) & saved_data_mat(:, 4) < bin_vec(bin_n + 1) );
    saved_means_bins(bin_n, :) = mean(saved_data_mat(curr_bini, :), 'omitnan');
    saved_ses_bins(bin_n, :) = std(saved_data_mat(curr_bini, :), [], 'omitnan')./sqrt(length(curr_bini));
end


%plotting
figure(1)
plot(saved_data_mat(:, 3), saved_data_mat(:, 1), 'O')
ylabel('peak response (dF/F)')
xlabel('session time (s)')
fig_wrapup(1, [])
ax = axis;
ax(4) = 0.5;
ax(3) = -0.2;
axis(ax);

figure(5)
errorbar(saved_means_bins(:, 3), saved_means_bins(:, 1), saved_ses_bins(:, 1))
ylabel('peak response (dF/F)')
xlabel('session time (s)')
fig_wrapup(5, [])
axis(ax);

figure(2)
plot(saved_data_mat(:, 3), saved_data_mat(:, 2), 'O')
ylabel('response area (dF/F)')
xlabel('session time (s)')
fig_wrapup(2, [])
axis(ax);

figure(6)
errorbar(saved_means_bins(:, 3), saved_means_bins(:, 2), saved_ses_bins(:, 2))
ylabel('response area (dF/F)')
xlabel('session time (s)')
fig_wrapup(6, [])
axis(ax);

figure(3)
plot(saved_data_mat(:, 4), saved_data_mat(:, 1), 'O')
ylabel('peak response (dF/F)')
xlabel('imaging time (s)')
fig_wrapup(3, [])
ax(2) = 2500;
axis(ax);

figure(7)
errorbar(saved_means_bins(:, 4), saved_means_bins(:, 1), saved_ses_bins(:, 1))
ylabel('peak response (dF/F)')
xlabel('imaging time (s)')
fig_wrapup(7, [])
axis(ax);


figure(4)
plot(saved_data_mat(:, 4), saved_data_mat(:, 2), 'O')
ylabel('response area (dF/F)')
xlabel('imaging time (s)')
fig_wrapup(4, [])
axis(ax);

figure(8)
errorbar(saved_means_bins(:, 4), saved_means_bins(:, 2), saved_ses_bins(:, 2))
ylabel('peak response (dF/F)')
xlabel('imaging time (s)')
fig_wrapup(8, [])
axis(ax);