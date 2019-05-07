clear all
close all

dataset_list_paths = [...
                      
                      {'C:\Data\Code\general_code_old\data_folder_lists\Janelia\MPPC_KC_set_final.xls'}...
                    ];
            
suppress_plots = 0;
[del, odor_names] = xlsread('C:\Data\Code\general_code_old\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList.xls', 1);
global color_vec;                
a = colormap('bone');
global greymap
greymap = flipud(a);
fly_n = 0;
script_name = mfilename;
resp_sizes_PMT = [];
resp_sizes_MPPC = [];
MPPC_colour = [0.9, 0.4, 0.6];
PMT_colour = [0.5, 0.6, 0.8];

for list_n = 1:size(dataset_list_paths, 1)
    curr_dir_list_path = dataset_list_paths{list_n, 1};
    [del, dir_list] = xlsread(curr_dir_list_path, 1);        %list of Suite2P results directories
    n_dirs = size(dir_list, 1);
    
    %loop to go through all experiment datasets listed in list file
    for dir_n = 1:n_dirs
               
        saved_an_results.scriptname = mfilename('fullpath');
        curr_dir = [dir_list{dir_n, 1}, '\'];
        curr_dir = manage_base_paths(curr_dir, 2);
       
        cd(curr_dir);
        tif_name = dir('*.tif');
        stack_obj = ScanImageTiffReader([curr_dir, tif_name(1).name]);
        [frame_time, zoom, n_chans, PMT_offsets] = SI_tif_info(stack_obj);
        
        tif_times = load([curr_dir, 'tif_time_stamps.mat']);           %reading in time stamps for each tif file recorded by raw_data_extracter
        tif_times = tif_times.time_stamps;
        [stim_mat, stim_mat_simple, column_heads, color_vec] = load_params_trains_modular(curr_dir, tif_times);    %reading in trial stimulus parameters after matching time stamps to F traces
               
        stim_frs = compute_stim_frs(stim_mat, 1, frame_time);
        stim_frs = stim_frs + round(.65./frame_time);
        
        %Reading in experimental parameters
        odor_list = unique(stim_mat_simple(:, 2) );
        n_odors = length(odor_list);
        odor_dur_list = unique(stim_mat_simple(:, 3) );
        n_od_durs = length(odor_dur_list);
        n_trains = max(stim_mat_simple(:, 11));
        
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
        [resp_sizes, sig_trace_mat, sig_trace_mat_old, sig_cell_mat, resp_areas] = cal_sig_responses_res(dff_data_mat, stim_mat, stim_mat_simple, curr_dir, frame_time);
        sig_cells = find(sum(sum(sig_cell_mat, 3), 2) > 0);         %list of all cells significant for any odor for any duration
        
        %separating PMT and MPPC response matrices
        n_frs_total = size(dff_data_mat, 1);
        dff_data_mat_PMT = dff_data_mat(1:2:n_frs_total, :, :);
        dff_data_mat_MPPC = dff_data_mat(2:2:n_frs_total, :, :);
        dff_data_mat_PMT_f = dff_data_mat_f(1:2:n_frs_total, :, :);
        dff_data_mat_MPPC_f = dff_data_mat_f(2:2:n_frs_total, :, :);
        raw_data_mat_PMT = raw_data_mat(1:2:n_frs_total, :, :);
        raw_data_mat_MPPC = raw_data_mat(2:2:n_frs_total, :, :);
        
        
        %normalising resp_matrices
        mean_PMT_sig = mean(raw_data_mat_PMT(1:100, 1, 1));
        dff_data_mat_PMT = dff_data_mat_PMT./mean_PMT_sig;
        mean_MPPC_sig = mean(raw_data_mat_MPPC(1:100, 1, 1));
        dff_data_mat_MPPC = dff_data_mat_MPPC./mean_MPPC_sig;
        
        %plotting response matrix differences
        diff_dff_data_mat = dff_data_mat_PMT - dff_data_mat_MPPC;
        
        %identifying sig responses
        for odor_n = 1:length(odor_list)
            odor_ni = odor_list(odor_n);
            curr_trs = find(stim_mat_simple(:, 2) == odor_ni);
            for cell_n = 1:n_cells
                ave_tr = mean(dff_data_mat_PMT(:, cell_n, curr_trs), 3, 'omitnan');
                ave_tr_f = movmean(ave_tr, 15, 'omitnan');
                max_val = max(ave_tr_f(stim_frs(1):(stim_frs(2) + round(2./frame_time)), 1), [], 1, 'omitnan');
                sd_baseline = std(ave_tr_f((stim_frs(1) - round(6./frame_time)):(stim_frs(1) - 1)), [], 'omitnan');
                
                if max_val./sd_baseline > 5
                    sig_cell_mat(cell_n, odor_ni) = 1;
                    curr_resp_sizes = squeeze(max(dff_data_mat_PMT((stim_frs(1)):(stim_frs(2) + round(2./frame_time)), cell_n, curr_trs)));
                    curr_resp_sizes_m = mean(curr_resp_sizes);
                    resp_sizes_PMT = [resp_sizes_PMT; curr_resp_sizes];
                    resp_sizes_MPPC = [resp_sizes_MPPC; squeeze(max(dff_data_mat_MPPC((stim_frs(1) ):(stim_frs(2) + round(2./frame_time)), cell_n, curr_trs)))];
                    
                    if suppress_plots == 0 %&& curr_resp_sizes_m > 1.5
                        figure(1)
                        %plotting traces
                        plot(squeeze(dff_data_mat_PMT((stim_frs(1) - round(5./frame_time)):(stim_frs(2) + round(5./frame_time)), cell_n, curr_trs(1)))...
                                                                                                                   , 'Color', PMT_colour, 'LineWidth', 0.5);
                        hold on
                        plot(squeeze(dff_data_mat_MPPC((stim_frs(1) - round(5./frame_time)):(stim_frs(2) + round(5./frame_time)), cell_n, curr_trs(1)) + 2.5)... 
                                                                                                                   , 'Color', MPPC_colour, 'LineWidth', 0.5);

                        set_xlabels_time(1, frame_time, 2)
                        ylabel('Fluorescence intensity (dF/F)')
                        ax_vals = axis;
                        ax_vals(3) = - 1;
                        axis(ax_vals);                                                                                       
                        hold off                                                                                       
                        fig_wrapup_nonums(1, script_name, 2); 
                        add_stim_bar(1, (stim_frs - (stim_frs(1) - round(5./frame_time))), [0.7, 0.7, 0.7])

                        figure(1)
                        keyboard
                        %del = input('press enter');
                        close figure 1
                    else
                    end
                    
                    
                else
                end
            end
            
        end
        
        figure(2)
        plot(resp_sizes_PMT, resp_sizes_MPPC, '.', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')
        hold on
        %fitting line
        p = polyfit(resp_sizes_PMT, resp_sizes_MPPC, 1);
        x1 = [0, 5];
        fit_vals = polyval(p, x1);
        plot(x1, fit_vals, '--', 'Color', [0.7, 0.7, 0.7]);
        hold off
        
        xlabel('PMT response size (dF/F)');
        ylabel('MPPC response size (dF/F)');
        axis([0, 5, 0, 5]);
        fig_wrapup(2, script_name);
        
        
        
        
        
        keyboard
    end
    
    
    
end