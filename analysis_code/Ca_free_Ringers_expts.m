clear all
close all

% global color_vec;                
color_vec = load('C:\Users\Mehrab\Google Drive\Backup\Stuff\CSHL\Glenn lab\Code\std_color_vec.txt');

direc_list_paths = [{'C:\Data\Data\Raw_data\dataset_lists\dataset_list_Cafree_Ringers_expt.xls'},...
                    {'C:\Data\Data\Raw_data\dataset_lists\dataset_list_Cafree_Ringers_expt - ctrl.xls'}];


for dir_list_n = 1:size(direc_list_paths, 2)
    direc_list_path = direc_list_paths{1, dir_list_n};
    [del, direc_list] = xlsread(direc_list_path, 1);
    n_direcs = size(direc_list, 1); 
    
    saved_resp_vecs = [];
    if isempty(findstr(direc_list_path, 'ctrl')) == 1
        if exist('C:\Data\Data\Analysed_data\Analysis_results\CaFreeRingersExpt\expt_traces.mat') == 2
            saved_resp_vecs = load('C:\Data\Data\Analysed_data\Analysis_results\CaFreeRingersExpt\expt_traces.mat');
            analyse = 0;
            saved_resp_vecs = saved_resp_vecs.saved_resp_vecs;
        else
            analyse = 1;
        end
    elseif isempty(findstr(direc_list_path, 'ctrl')) == 0
        if exist('C:\Data\Data\Analysed_data\Analysis_results\CaFreeRingersExpt\ctrl_traces.mat') == 2
            saved_resp_vecs = load('C:\Data\Data\Analysed_data\Analysis_results\CaFreeRingersExpt\ctrl_traces.mat');
            analyse = 0;
            saved_resp_vecs = saved_resp_vecs.saved_resp_vecs;
        else
            analyse = 1;
        end
    else
    end


    if analyse == 1
        for direc_n = 1:n_direcs
            curr_dir = direc_list{direc_n, 1};
            curr_dir = [curr_dir, '\'];
            remove_small_tifs(curr_dir);
            dir_contents = dir_date_sorted(curr_dir, '*.tif');
            [stim_mat, stim_mat_simple, column_heads] = load_params_trains(curr_dir, []);

            n_trials = size(stim_mat_simple, 1);
            mean_resps = zeros(n_trials, 1) + nan;

            for trial_n = 1:n_trials
                stack_obj = ScanImageTiffReader([curr_dir, dir_contents(trial_n).name]);
                stack = stack_obj.data();
                stack = permute(stack,[2 1 3]);
                n_frames = size(stack, 3);
                [frame_time, zoom, n_chans] = SI_tif_info(stack_obj);
                stim_fr = floor(stim_mat_simple(trial_n, 7)./frame_time);
                stim_end_fr = floor( (stim_mat_simple(trial_n, 7) + stim_mat_simple(trial_n, 3))./frame_time);

                if trial_n == 1
                    [del, sig_resp_frame, del2, smple_fr] = highlight_resp_pix(1, stack, [stim_fr, stim_end_fr], 1, frame_time, 0);
                elseif trial_n > 1    
                    [del, sig_resp_frame, del2, smple_fr] = highlight_resp_pix(1, stack, [stim_fr, stim_end_fr], 1, frame_time, 1);
                else
                end
                sig_resp_frame(sig_resp_frame == 0) = nan;                              %excluding non-responsive pixels from calculation
                %mean_resps(trial_n, 1) = sum(sum(sig_resp_frame, 'omitnan'), 'omitnan');
                mean_resps(trial_n, 1) = length(find(isnan(sig_resp_frame) == 0))./sum(sum(smple_fr));
                disp(['done with trial ' int2str(trial_n) ' of ' int2str(n_trials)])        
            end
            figure(2)
            odor_list = unique(stim_mat_simple(:, 2));
            for odor_n = 1
                odor_ni = odor_list(odor_n);
                curr_trs = find(stim_mat_simple(:, 2) == odor_ni);

                plot(mean_resps(curr_trs), 'Color', color_vec(odor_n, :))
                hold on
                saved_resp_vecs{direc_n, odor_n} = mean_resps(curr_trs);
                ylabel('frac. responsive pixels');
            end

            if isempty(findstr(direc_list_path, 'ctrl')) == 1
                save(['C:\Data\Data\Analysed_data\Analysis_results\CaFreeRingersExpt\expt_traces.mat'], 'saved_resp_vecs');
            elseif isempty(findstr(direc_list_path, 'ctrl')) == 0
                save(['C:\Data\Data\Analysed_data\Analysis_results\CaFreeRingersExpt\ctrl_traces.mat'], 'saved_resp_vecs');
            else 
            end
        end
    else
    end
end

%loading in data and plotting
saved_resp_vecs = load('C:\Data\Data\Analysed_data\Analysis_results\CaFreeRingersExpt\expt_traces.mat');
saved_resp_vecs = saved_resp_vecs.saved_resp_vecs;
saved_resp_vecs_ctrl = load('C:\Data\Data\Analysed_data\Analysis_results\CaFreeRingersExpt\ctrl_traces.mat');
saved_resp_vecs_ctrl = saved_resp_vecs_ctrl.saved_resp_vecs;

transition_trial_ns = [6, 36, 41];
trace_mat = [];
for expt_n = 1:2:3
    curr_trace = saved_resp_vecs{expt_n, 1};
    curr_trace = curr_trace(1:45);
    trace_mat = [trace_mat, curr_trace];
end
mean_trace = mean(trace_mat, 2, 'omitnan');
sd_trace = std(trace_mat, 0, 2, 'omitnan');

figure(3)
shadedErrorBar([(1:45).*2], mean_trace, sd_trace, {'Color', [0.3, 0.7, 0.3]});
hold on
shadedErrorBar([(6:40).*2], mean_trace(6:40), sd_trace(6:40), {'Color', [0.7, 0.3, 0.3]});

ylabel('frac. responsive pixels');
xlabel('time (min)')

%control traces
trace_mat = [];
for expt_n = 1:2
    curr_trace = saved_resp_vecs_ctrl{expt_n, 1};
    curr_trace = curr_trace(1:end);
    trace_mat = pad_n_concatenate(trace_mat, curr_trace, 2, nan);
end

mean_trace = mean(trace_mat, 2, 'omitnan');
sd_trace = std(trace_mat, 0, 2, 'omitnan');

figure(3)
%shadedErrorBar([(1:40).*2], mean_trace(1:40), sd_trace(1:40), {'Color', [0.6, 0.6, 0.6]}, 0);
fig_wrapup(3)
x = [0.265 0.265];
y = [0.95 0.78];
a = annotation('textarrow',x,y,'String','Ca-free Ringer''s wash-in');

x = [0.76 0.76];
y = [0.46 0.29];
b = annotation('textarrow',x,y, 'String', 'Ca-Ringer''s wash-in');

x = [0.7 0.7];
y = [0.5 0.33];
b = annotation('textarrow',x,y, 'String', 'removed glial sheath');
