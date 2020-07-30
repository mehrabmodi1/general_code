clear all
close all


path_base = 'C:\Data\Data\Analysed_data\data_sharing\param_fitting\';
frame_time = 0.0999;    %in s

cd(path_base);
dir_contents = dir;
dir_contents(1:2) = [];
cell_type_list = [];
fly_n_list = [];
for type_n = 1:size(dir_contents, 1)
    cell_type_list = [cell_type_list; {dir_contents(type_n).name}];
    cd([path_base, cell_type_list{type_n}]);
    n_flies = dir();
    fly_n_list = [fly_n_list; (size(n_flies, 1) - 2)];
end

a = colormap('bone');
greymap = flipud(a);
model_traces_all = [];
resp_traces_all = [];
cell_params_all = [];

color_vecs = cbrewer('qual', 'Dark2' , size(cell_type_list, 1), 'cubic');
type_name_vec = [];
params_saved_all = [];
for cell_type_n = 1:size(cell_type_list, 1)
    cell_type_name = cell_type_list{cell_type_n};
    cell_type_namei = findstr(cell_type_name, 'KC');
    cell_type_name = cell_type_name(1:(cell_type_namei - 1));
    type_name_vec = [type_name_vec; {cell_type_name}];
    
    n_flies = fly_n_list(cell_type_n, 1);
    
    curr_params_saved = [];
    for fly_n = 1:n_flies
        an_path = [path_base, cell_type_list{cell_type_n}, '\fly', num2str(fly_n), '\'];

        %reading in KC response params fitted with Herve's program
        fit_params = readNPY([an_path, 'fit_params.npy']);

        %reading in KC response traces used by python script to generate fits along with other metadata
        fit_tr_list = load([an_path, 'tr_list.mat']);       %this is the list of trials with only long dur, single od pulse. Used for the KC param fits
        fit_tr_list = fit_tr_list.tr_list;
        traces_frompy = readNPY([an_path, 'fit_data.npy']);
        dff_data_mat = load([an_path, 'dFF_data.mat']);
        dff_data_mat = dff_data_mat.dff_data_mat_f;
        stim_mat = load([an_path, 'stim_mat.mat']);
        stim_mat = stim_mat.stim_mat;

        [stim_mat_simple, column_heads] = generate_stim_mat_simple(stim_mat);
        dur_olf1_col_n = find_stim_mat_simple_col('duration', column_heads);        %identifying relevant column number in stim_mat_simple
        dur_olf2_col_n = find_stim_mat_simple_col('duration_olf2', column_heads);   %identifying relevant column number in stim_mat_simple
        olf1_od_col_n = find_stim_mat_simple_col('odor_n', column_heads);        %identifying relevant column number in stim_mat_simple
        olf2_od_col_n = find_stim_mat_simple_col('odour_olf2', column_heads);        %identifying relevant column number in stim_mat_simple
        od_col_ns = [olf1_od_col_n, olf2_od_col_n];
        dur_col_ns = [dur_olf1_col_n, dur_olf2_col_n];

        [resp_sizes, sig_trace_mat, sig_trace_mat_old, sig_cell_mat, resp_areaundercurves] = cal_sig_responses_res_modular(dff_data_mat, stim_mat, stim_mat_simple, frame_time, od_col_ns, dur_col_ns);

        od_list_olf1 = unique(stim_mat_simple(:, olf1_od_col_n));


        %computing fits
        for od_n = 1:length(od_list_olf1)
            od_ni = od_list_olf1(od_n);
            curr_trs = find(stim_mat_simple(:, olf1_od_col_n) == od_ni);      %list of trials with current 
            param_od_ns = [od_n, nan];      %odor indices for olf1 and olf2; ie not odor numbers, but their positions in odor lists.     
            try
                [model_traces, resp_traces, cell_params, stim_frs] = get_fit_traces(curr_trs, stim_mat, od_list_olf1, fit_params, frame_time, dff_data_mat, sig_cell_mat);
            catch
                [model_traces, resp_traces, cell_params] = get_fit_traces(curr_trs, stim_mat, od_list_olf1, fit_params, frame_time, dff_data_mat, sig_cell_mat);
            end


            if isempty(resp_traces) == 0
                model_traces_all = pad_n_concatenate(model_traces_all, model_traces, 2, nan);
                resp_traces_all = pad_n_concatenate(resp_traces_all, resp_traces, 2, nan);
                cell_params_all = pad_n_concatenate(cell_params_all, cell_params, 2, nan);
            else
            end
        end   

        
        %re-shaping and saving fit_params to a single matrix, keeping track of
        %cell number and cell type as params 9 and 10

        %adding on params 9 (cell number) and 10 (dataset list number)
        if fly_n == 1
            start_cell_n = 1;
        elseif fly_n > 1
            start_cell_n = end_cell_n + 1;
        else
        end
        end_cell_n = start_cell_n + size(fit_params, 2) - 1;
        nine_mat = repmat( (start_cell_n:end_cell_n), size(fit_params, 1), 1, 1 );
        fit_params(:, :, 9) = nine_mat;
        fit_params(:, :, 10) = cell_type_n;

        %re-arranging fit_params so that each odor-cell pair's parameter set is
        %a separate row and then accumulating across flies
        for od_n = 1:size(fit_params, 1)
            curr_params = [squeeze(fit_params(od_n, :, :)), repmat(od_n, size(fit_params, 2), 1)];
            curr_params_saved = [curr_params_saved; curr_params];
        end
             
    end
     
    %sorting traces by centers of mass of model traces
    centers_of_mass = zeros(size(model_traces_all, 2), 1) + nan;
    for cell_n = 1:size(model_traces_all, 2)
        cOm = centerOfMass(model_traces_all(:, cell_n));
        centers_of_mass(cell_n, 1) = cOm(1, 1);
    end
    model_traces_all = [centers_of_mass'; model_traces_all];
    model_traces_all = sortrows(model_traces_all')';
    model_traces_all(1, :) = [];
    resp_traces_all = [centers_of_mass'; resp_traces_all];
    resp_traces_all = sortrows(resp_traces_all')';
    resp_traces_all(1, :) = [];


    %plotting raw data
    n_frs = size(model_traces_all, 1) - 10;
    fig_h1 = figure('Name', 'dF/F response traces');        
    imagesc(resp_traces_all(1:n_frs, :)')
    colormap(greymap)
    set_xlabels_time(fig_h1, frame_time, 20);
    ylabel('cell-odor pairs');
    fig_wrapup(fig_h1, []);
    add_stim_bar(fig_h1, stim_frs, [0.65, 0.65, 0.65]);

    fig_h2 = figure('Name', 'modelled response traces');        
    imagesc(model_traces_all')
    colormap(greymap)
    set_xlabels_time(fig_h2, frame_time, 20);
    ylabel('cell-odor pairs');
    fig_wrapup(fig_h2, []);
    add_stim_bar(fig_h2, stim_frs, [0.65, 0.65, 0.65]);
    
    %accumulating params across cell-types
    params_saved_all = [params_saved_all; curr_params_saved];
end

param_name_list = [{'Aon'}, {'Asus'}, {'Aoff'}, {'T1on'}, {'T1off'}, {'T2on'}, {'T2off'}];

%1. plotting distributions of each parameter for each cell-type
param_max_vals = max(params_saved_all, [], 'omitnan');
param_max_vals = param_max_vals(2:8);
param_max_vals = [param_max_vals; [3, 3, 3, 5, 5, 5, 5]];
param_max_vals = min(param_max_vals, [], 'omitnan');
keyboard
close all
n_bins = 10;
for cell_type_n = 1:size(cell_type_list, 1)
    curr_cells = find(params_saved_all(:, 10) == cell_type_n);
    curr_param_mat = params_saved_all(curr_cells, :);
    
    
    for param_n = 2:8
        bin_vec = 0:(param_max_vals./n_bins):param_max_vals((param_n - 1));
        curr_par_vals = curr_param_mat(:, param_n);
        counts = histcounts(curr_par_vals, bin_vec);
        counts = counts./sum(counts);
        figure( (param_n - 1) )
        plot(bin_vec(2:end), counts, 'Color', color_vecs(cell_type_n, : ), 'lineWidth', 2.5, 'DisplayName', type_name_vec{cell_type_n});
        if param_n < 5
            xlabel([param_name_list{param_n - 1}, ' values (dF/F)'])
        elseif param_n >= 5
            xlabel([param_name_list{param_n - 1}, ' values (s)'])
        else
        end
        ylabel('fraction of all counts')
        fig_wrapup((param_n - 1), []);
        hold on
        
        if cell_type_n == size(cell_type_list, 1)
            legend('show');
        else
        end
        
    end
end


%clustering cell-od pairs based on corrcoefs between parameter vecs
data_mat = params_saved_all(:, 5:8);    %first amplit parameter is for pre-stimulus baseline and is always close to 0
label_mat = params_saved_all(:, [1:4, 9:11]);  %matrix of cell number, cell type number and odor number labels

%getting rid of all-0 rows
row_sum = sum(data_mat, 2);
del = find(row_sum == 0);
data_mat(del, :) = [];
label_mat(del, :) = [];

n_clusts = 4;

[data_mat_sorted, data_mat_sorted_n, label_mat_sorted, clust_ids_sorted, link_map] = cluster_euc_hier(data_mat, label_mat, n_clusts);

corr_mat = corrcoef(data_mat_sorted_n', 'Rows', 'complete');

figure(8)
imagesc(data_mat_sorted)

figure(9)
imagesc(corr_mat)

figure(10)
imagesc(label_mat_sorted, [0, 3])


%plotting modelled response trace for each cluster
cell_type_mat = zeros(n_clusts, size(cell_type_list, 1));
clust_traces = [];
for clust_n = 1:n_clusts
    curr_train = [0, 30];
    curr_members = find(clust_ids_sorted == clust_n);
    
    %keeping track of fractions of cell types in current cluster
    curr_cell_types = label_mat(curr_members, 2);
    for cell_type_n = 1:size(cell_type_list, 1)
        cell_type_mat(clust_n, cell_type_n) = length(find(curr_cell_types == cell_type_n));
    end
    
    curr_params_taus = mean(data_mat_sorted(curr_members, :), 1);
    curr_params_amplits = mean(label_mat_sorted(curr_members, 2:4), 1);
    curr_params = [0, curr_params_amplits, curr_params_taus];
    curr_params = mean(curr_params, 1);
        
    %getting model trace for average cell of current cluster
    [model_trace, stim_frs] = fit_KC_response2(curr_train, curr_params, frame_time, 1400, 10);
    
    clust_traces = [clust_traces, model_trace];
   
    
end
figure(11)
plot(clust_traces)


%PICK UP THREAD HERE
%1. get params for the same stim durations across cell types and compute
%averaged cluster response trace for cluster types
%2. get larger numbers of repeats, extract params for different sets of
%repeats and compare trace reconstructions with un-fitted response traces
%3. statistical testing comparing different KC subtypes


keyboard

function [model_traces, resp_traces, cell_params, stim_frs] = get_fit_traces(curr_trs, stim_mat, od_list_olf1, fit_params, frame_time, dff_data_mat, sig_cell_mat)
%determining if this is an olf1 or an olf2 odor
curr_od = stim_mat(curr_trs(1)).odours;
if curr_od <= 12    %case where it is an olf1 odor
    curr_train = stim_mat(curr_trs(1)).pulse_train;
elseif curr_od > 12 %case where it is an olf2 odor
    curr_train = stim_mat(curr_trs(1)).pulse_train_olf2;
else
end

n_frames = size(dff_data_mat, 1) - sum(isnan(dff_data_mat(:, 1, curr_trs(1)))); 
stim_latency = stim_mat(curr_trs(1)).stimLatency;

od_n = find(od_list_olf1 == curr_od);   %identifying position in odorlist for current odor to correctly choose fit params
curr_sig_cells = sig_cell_mat(:, curr_od);
model_traces = [];
resp_traces = [];
cell_params = [];
for cell_n = 1:size(fit_params, 2)
    curr_resp_trace = mean(dff_data_mat(:, cell_n, curr_trs), 3, 'omitnan');
    curr_cell_params = squeeze(fit_params(od_n, cell_n, :));
    %checking if current cell has a significant response trial in curr_trs
    %and skipping if not
    if curr_sig_cells(cell_n) < 1
        continue
    elseif isnan(curr_resp_trace(1)) == 1
        continue
    elseif sum(curr_cell_params) == 0
        continue
    else
    end
    
    [model_trace, stim_frs] = fit_KC_response(curr_train, curr_cell_params, frame_time, n_frames, stim_latency);
    
    %logging data
    model_traces = [model_traces, model_trace];
    resp_traces = [resp_traces, curr_resp_trace];
    cell_params = [cell_params, curr_cell_params;];
    
end


function [model_trace, stim_frs] = fit_KC_response(curr_train, curr_cell_params, frame_time, n_frames, stim_latency)
a0 = curr_cell_params(1);
a1 = curr_cell_params(2);
a2 = curr_cell_params(3) ;
a3 = curr_cell_params(4);
t0 = curr_cell_params(5);
t1 = curr_cell_params(6);
t2 = curr_cell_params(7);
t3 = curr_cell_params(8);

n_pulses = size(curr_train, 1);

t = ((1:1:n_frames) * frame_time)';
model_trace = zeros(n_frames, 1);
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
    stim_frs = (curr_train + stim_latency)./frame_time;
%         plot_fit_trace(1, curr_resp_trace, curr_model_trace, stim_frs, frame_time);
%         keyboard
%         close figure 1


end
    
end  
end


function [model_trace, stim_frs] = fit_KC_response2(curr_train, curr_cell_params, frame_time, n_frames, stim_latency)
    a0 = curr_cell_params(1);
    a1 = curr_cell_params(2);
    a2 = curr_cell_params(3) ;
    a3 = curr_cell_params(4);
    t0 = curr_cell_params(5);
    t1 = curr_cell_params(6);
    t2 = curr_cell_params(7);
    t3 = curr_cell_params(8);

    n_pulses = size(curr_train, 1);
    
    t = ((1:1:n_frames) * frame_time)';
    model_trace = zeros(n_frames, 1);
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
        stim_frs = (curr_train + stim_latency)./frame_time;
%         plot_fit_trace(1, curr_resp_trace, curr_model_trace, stim_frs, frame_time);
%         keyboard
%         close figure 1
      
        
    end
    
end

function [] = plot_fit_trace(fig_n, resp_trace, model_trace, stim_frs, frame_time)

%plotting response with fit overlaid
figure(fig_n)
plot(resp_trace);
hold on
plot(model_trace', 'r')
set_xlabels_time(1, frame_time, 20);
ylabel('response amplitude (dF/F)')
fig_wrapup(fig_n, []);
add_stim_bar(1, stim_frs, [0.65, 0.65, 0.65]);
end
