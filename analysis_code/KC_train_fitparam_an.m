clear all
close all

path_base = 'C:\Data\Data\Analysed_data\data_sharing\KC_long_trace\';

n_flies = 4;
frame_time = 0.0999;    %in s

saved_model_traces = [];
saved_resp_traces = [];
%cell_params_all = [];
for fly_n = 1:n_flies
    an_path = [path_base, 'fly', num2str(fly_n), '\'];
    
    %reading in KC response params fitted with Herve's program
    fit_params = readNPY([an_path, 'fit_params.npy']);
    
    %reading in KC response traces used by python script to generate fits along with other metadata
    fit_tr_list = load([an_path, 'tr_list.mat']);       %this is the list of trials with only long dur, single od pulse. Used for the KC param fits
    fit_tr_list = fit_tr_list.tr_list;
    traces_frompy = readNPY([an_path, 'fit_data.npy']);
    dff_data_mat_trimmed = load([an_path, 'dFF_data.mat']);
    dff_data_mat_trimmed = dff_data_mat_trimmed.dff_data_mat_f;
    stim_mat_trimmed = load([an_path, 'stim_mat.mat']);
    stim_mat_trimmed = stim_mat_trimmed.stim_mat;
    dff_data_mat = load([an_path, 'dFF_data_mat_orig.mat']);
    dff_data_mat = dff_data_mat.dff_data_mat_f_orig;
    stim_mat = load([an_path, 'stim_mat_orig.mat']);
    stim_mat = stim_mat.stim_mat_orig;
    
    [stim_mat_simple, column_heads] = generate_stim_mat_simple(stim_mat);
    dur_olf1_col_n = find_stim_mat_simple_col('duration', column_heads);        %identifying relevant column number in stim_mat_simple
    dur_olf2_col_n = find_stim_mat_simple_col('duration_olf2', column_heads);   %identifying relevant column number in stim_mat_simple
    olf1_od_col_n = find_stim_mat_simple_col('odor_n', column_heads);        %identifying relevant column number in stim_mat_simple
    olf2_od_col_n = find_stim_mat_simple_col('odour_olf2', column_heads);        %identifying relevant column number in stim_mat_simple
    od_col_ns = [olf1_od_col_n, olf2_od_col_n];
    dur_col_ns = [dur_olf1_col_n, dur_olf2_col_n];
    
    [resp_sizes, sig_trace_mat, sig_trace_mat_old, sig_cell_mat, resp_areaundercurves] = cal_sig_responses_res_modular(dff_data_mat, stim_mat, stim_mat_simple, frame_time, od_col_ns, dur_col_ns);
    
    [stim_mat_simple_trimmed, column_heads] = generate_stim_mat_simple(stim_mat_trimmed);
    od_list_trimmed = unique(stim_mat_simple_trimmed(:, olf1_od_col_n));
    od_list_olf1 = od_list_trimmed(od_list_trimmed <= 12);          %list of olf1 odors for which we have fitted params
    od_list_olf2 = od_list_trimmed(od_list_trimmed >= 12) - 12;     %list of olf2 odors for which we have fitted params
    stim_mat_simple(isnan(stim_mat_simple)) = 0;    
    dur_list_olf1 = unique(stim_mat_simple(:, dur_col_ns(1)));
    dur_list_olf1 = dur_list_olf1(dur_list_olf1 > 0.1);
    dur_list_olf2 = unique(stim_mat_simple(:, dur_col_ns(2)));
    dur_list_olf2 = dur_list_olf2(dur_list_olf2 > 0.1);
    
    %defining a lookup table with three columns to allow referencing the fit_params matrix given olf1 and olf2 od numbers
    param_mat_od_list = unique(stim_mat_simple_trimmed(:, od_col_ns(1)));
    od_n_lookup = zeros(length(param_mat_od_list), 3);      %initialising lookup table
    od_n_lookup(:, 3) = 1:1:length(param_mat_od_list);
    od_n_lookup(:, 1:2) = [param_mat_od_list, param_mat_od_list];
    del1 = find(od_n_lookup(:, 1) > 12);        %indices of olf1_ods
    od_n_lookup(del1, 1) = nan;
    del2 = find(od_n_lookup(:, 2) <= 12);       %indices of olf2_ods
    od_n_lookup(del2, 2) = nan;
    od_n_lookup(del1, 2) = od_n_lookup(del1, 2) - 12;
    
    %computing fits for handover trials and rand train trials, if any
    for odor_n_olf1 = 1:length(od_list_olf1)
        olf1_od_ni = od_list_olf1(odor_n_olf1);
        %checking if current olf1_od has params fitted for it and skipping if not
        del = find(od_n_lookup(:, 1) == olf1_od_ni);
        if isempty(del) == 1
            continue
        else
            olf1_param_mat_od_n = od_n_lookup(del, 3);
        end
                
        
        for odor_n_olf2 = 1:length(od_list_olf2)
            olf2_od_ni = od_list_olf2(odor_n_olf2);
            
            %checking if current olf2_od has params fitted for it and skipping if not
            del = find(od_n_lookup(:, 2) == olf2_od_ni);
            if isempty(del) == 1
                continue
            else
                olf2_param_mat_od_n = od_n_lookup(del, 3);
            end

            for dur_n_olf1 = 1:length(dur_list_olf1)
                olf1_dur = dur_list_olf1(dur_n_olf1);
                for dur_n_olf2 = 1:length(dur_list_olf2)
                    olf2_dur = dur_list_olf2(dur_n_olf2);
                    curr_trs = find(stim_mat_simple(:, od_col_ns(1)) == olf1_od_ni & stim_mat_simple(:, od_col_ns(2)) == olf2_od_ni &...
                                        stim_mat_simple(:, dur_col_ns(1)) == olf1_dur & stim_mat_simple(:, dur_col_ns(2)) == olf2_dur);
                    
                    if isempty(curr_trs) == 1
                        continue
                    else
                    end            
                                    
                    olf1_train = stim_mat(curr_trs(1)).pulse_train;
                    olf2_train = stim_mat(curr_trs(1)).pulse_train_olf2;
                    olf2_train = olf2_train + stim_mat(curr_trs(1)).rel_stimLatency_olf2;
                    olf_trains = [{olf1_train}; {olf2_train}];
                    
                    %making sure current trials are handover or odor train trials
                    if isnan(olf2_train) == 1 | max(olf1_train) == 0.1
                        %case when only one odor was delivered in curr_trs
                        keyboard
                        if max([size(olf1_train, 1), size(olf2_train, 1)]) > 1
                            %multi-pulse train on olf1 or olf2
                        else
                            %case where only a single pulse of a single odor was delivered on curr_trs
                            continue
                        end
                    else
                    end
                    
                    
                    n_frames = size(dff_data_mat, 1) - sum(isnan(dff_data_mat(:, 1, curr_trs(1)))); 
                    curr_resp_traces = squeeze(mean(dff_data_mat(1:n_frames, :, curr_trs), 3, 'omitnan'));
                    stim_latency = stim_mat(curr_trs(1)).stimLatency;
                    
                    for cell_n = 1:size(dff_data_mat, 2)
                        od_n_olf2i = olf2_od_ni + 12;
                        
                        %making sure current cell is a sig responder for at least one of the current odors
                        if sig_cell_mat(cell_n, olf1_od_ni) == 0 & sig_cell_mat(cell_n, od_n_olf2i) == 0
                            continue
                        else
                        end
                        
                        curr_resp_trace = curr_resp_traces(:, cell_n);
                        curr_cell_params = zeros(2, 8);
                        curr_cell_params(1, :) = squeeze(fit_params(olf1_param_mat_od_n, cell_n, :));
                        curr_cell_params(2, :) = squeeze(fit_params(olf2_param_mat_od_n, cell_n, :));
                       
                        %skipping cells with no fit params (all zeros) for either odor
                        if sum(curr_cell_params(1, :)) == 0 | sum(curr_cell_params(2, :)) == 0
                            continue
                        else
                        end
                        %calling fitting function with olf1_odn olf2_od_n, olf1_param_odn, and olf2_param_odn
                        [model_trace, stim_frs] = fit_KC_response(olf_trains, curr_cell_params, frame_time, n_frames, stim_latency, curr_resp_trace);
                        saved_model_traces = pad_n_concatenate(saved_model_traces, model_trace, 2, nan);
                        
                        saved_resp_traces = pad_n_concatenate(saved_resp_traces, curr_resp_trace, 2, nan);
                        
                    end
                    
                end
            end
        end
    end
    
    
        
end




function [model_trace, stim_frs] = fit_KC_response(olf_trains, curr_cell_params, frame_time, n_frames, stim_latency, curr_resp_trace)
    

    model_trace = zeros(n_frames, 1);
    for olf_train_n = 1:2
        curr_train = olf_trains{olf_train_n};
        n_pulses = size(curr_train, 1);
        t = ((1:1:n_frames) * frame_time)';
        
        %a0 = curr_cell_params(olf_train_n, 1);
        a0 = 0;
        a1 = curr_cell_params(olf_train_n, 2);
        a2 = curr_cell_params(olf_train_n, 3) ;
        a3 = curr_cell_params(olf_train_n, 4);
        t0 = curr_cell_params(olf_train_n, 5);
        t1 = curr_cell_params(olf_train_n, 6);
        t2 = curr_cell_params(olf_train_n, 7);
        t3 = curr_cell_params(olf_train_n, 8);
        
        for pulse_n = 1:n_pulses
            curr_model_trace = zeros(n_frames, 1);

            curr_pulse = curr_train(pulse_n, :);
            curr_train_time = sum(sum(curr_train(1:(pulse_n - 1), :)));
            
            on_time = stim_latency + curr_train_time + curr_pulse(1, 1);                       %time from beginning of train to onset of current pulse
            off_time = on_time + curr_pulse(1, 2);   %time from beginning of train to off of current pulse

            %boolean vectors that define the three epochs of the responses
            %to each pulse
            domain0 = t <= on_time;
            domain1 = (t > on_time) & (t <= off_time);
            domain2 = t > off_time;


            %constructing model response trace generated by current pulse
            curr_model_trace = curr_model_trace + domain0.*(a0);
            curr_model_trace = curr_model_trace + domain1.*(a1 .* (1.0 - exp(-(t - on_time) / t0))...
                                                            + a2 .* (1.0 - exp(-(t - on_time) / t1))...
                                                            + a0);
            %computing a scalar adjustment factor                                                        
            end1 = (a1 * (1.0 - exp(-on_time / t0))...
                    + a2 * (1.0 - exp(-on_time / t1))...
                    + a0);                       
            %adding the odor off part of the modelled response
            
            curr_model_trace = curr_model_trace + domain2.*(a3 .* exp(-(t - off_time) / t2)...
                                                            + (end1 - a3) .* exp(-(t - off_time) / t3));  
            
            model_trace = model_trace + curr_model_trace;       %adding up model responses for multiple odor pulses assuming linear summation
            
        end
        
        stim_frs{1} = (olf_trains{1} + stim_latency)./frame_time;
        stim_frs{2} = (olf_trains{2} + stim_latency)./frame_time;
        
        
    end
        
    plot_fit_trace(1, curr_resp_trace, model_trace, stim_frs, frame_time);
    keyboard
    close figure 1
    
end


function [] = plot_fit_trace(fig_n, resp_trace, model_trace, stim_frs, frame_time)

%plotting response with fit overlaid
figure(fig_n)
plot(resp_trace);
hold on
plot(model_trace', 'r')
set_xlabels_time(1, frame_time, 10);
fig_wrapup(fig_n, []);
add_stim_bar_modular(1, stim_frs, [0.65, 0.65, 0.65]);
end
