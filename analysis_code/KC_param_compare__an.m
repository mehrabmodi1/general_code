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

for cell_type_n = 1:size(cell_type_list, 1)
    n_flies = fly_n_list(cell_type_n, 1);
    
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
    keyboard
end

%plotting fit param distribution for off decay time or t3
fig_h3 = figure('Name', 'distribution of decay taus3');      
t3_vec = cell_params_all(8, :);         %multiplying by frame time to convert from n_frames to s
[counts, centers] = hist(t3_vec, 8);
counts = counts./sum(counts);
keyboard
bar(centers, counts)
xlabel('tau2 (s)')
ylabel('fraction of counts')
fig_wrapup(fig_h3, []);

%plotting fit param distribution for off decay time or t2
fig_h4 = figure('Name', 'distribution of decay taus2');      
t2_vec = cell_params_all(7, :);         %multiplying by frame time to convert from n_frames to s
[counts, centers] = hist(t2_vec, 8);
counts = counts./sum(counts);
bar(centers, counts)
xlabel('tau2 (s)')
ylabel('fraction of counts')
fig_wrapup(fig_h4, []);

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
    
    [model_trace, stim_frs] = fit_KC_response(curr_train, curr_cell_params, frame_time, n_frames, stim_latency, curr_resp_trace);
    
    %logging data
    model_traces = [model_traces, model_trace];
    resp_traces = [resp_traces, curr_resp_trace];
    cell_params = [cell_params, curr_cell_params;];
    
end

function [model_trace, stim_frs] = fit_KC_response(curr_train, curr_cell_params, frame_time, n_frames, stim_latency, curr_resp_trace)
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
