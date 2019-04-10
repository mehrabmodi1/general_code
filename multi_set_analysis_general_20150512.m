clear all
close all
list_direc = ['C:\Data\CSHL\dataset_list_Toshi_KC_DN_Led_20150708.txt'];
fid = fopen(list_direc);


saved_cells_b = [];
saved_cells_t = [];

saved_cells_b_shock = [];
saved_cells_t_shock = [];

%loop to go through all experiment datasets listed in list file
while 1
    
    direc = fgetl(fid);
       
    if ischar(direc) ~= 1
        break
    else
    end
    
    direc = [direc '\'];
    %checking if baseline and treatment trials were split into two sets
    if isdir([direc 'split_set']) == 1
        split_set = 1;
    else
        split_set = 0;
    end
    
    
    if split_set == 0
        dataset = load([direc 'expt.mat']);
        dataset = dataset.data;
        raw_data_mat = load([direc 'expt_raw_traces.mat']);
        raw_data_mat = raw_data_mat.raw_data_mat;
    elseif split_set == 1
        dataset = load([direc 'expt.mat']);
        dataset = dataset.data;
        baseline_data_mat = load([direc 'split_set\expt_raw_traces_baseline.mat']);
        baseline_data_mat = baseline_data_mat.raw_data_mat_baseline;
        treatment_data_mat = load([direc 'split_set\expt_raw_traces_treatment.mat']);
        treatment_data_mat = treatment_data_mat.raw_data_mat_treatment;
        
        %concatenating the two datasets for analysis
        n_frames = size(baseline_data_mat, 1);
        n_cells_b = size(baseline_data_mat, 2);
        n_cells_t = size(treatment_data_mat, 2);
        n_trials_b = size(baseline_data_mat, 3);
        n_trials_t = size(treatment_data_mat, 3);
        n_cells_tot = n_cells_b + n_cells_t;
        n_trials_tot = n_trials_b + n_trials_t;
        raw_data_mat = zeros(n_frames, n_cells_tot, n_trials_tot) + nan;
        
        raw_data_mat(1:n_frames, 1:n_cells_b, 1:n_trials_b) = baseline_data_mat;
        raw_data_mat(1:n_frames, (n_cells_b + 1):(n_cells_b + n_cells_t), (n_trials_b + 1):(n_trials_b + n_trials_t)) = treatment_data_mat;
        
    end
    
    %calculating dF/F traces and creating the sparse, 4-D, nan-filled
    %dff_data_mat 
    [dff_data_mat, stim_mat] = cal_dff_traces(raw_data_mat, dataset);
    clear raw_data_mat
    
    %identifying sig responses on a single trial basis, and then sig
    %responder cells
    [resp_areas, sig_trace_mat, sig_cell_mat] = cal_sig_responses(dataset, dff_data_mat, stim_mat);
    
    
    
    %PLOTTING
    n_odors = max(stim_mat(:, 1));
    
    %plotting trial-averaged response traces
    figure(1)
    for odor_n = 1:n_odors
        odor_name = dataset(1).stim.odourNames(odor_n).odour;
        subplot(ceil(n_odors/2), ceil(n_odors/2), odor_n)
        imagesc(squeeze(nanmean(dff_data_mat(:, :, :, odor_n), 3))', [0, 3]);
        title(['Trial-averaged traces for ' odor_name])
        xlabel('frame number')
        ylabel('cell number')
        
    end
    set(gcf, 'Color', 'w')
    
    %Plotting response areas over all trials
    figure(2)
    for odor_n = 1:n_odors
        odor_trs = find(stim_mat(:, 1) == odor_n);
        
        odor_name = dataset(1).stim.odourNames(odor_n).odour;
        subplot(ceil(n_odors/2), ceil(n_odors/2), odor_n)
        imagesc(resp_areas(:, odor_trs), [0, 30])
        
        title(['Area under the peak for ' odor_name])
        xlabel('trial number')
        ylabel('cell number')
        
    end
    set(gcf, 'Color', 'w')
    
    %Plotting response areas over trials, with and without elec stim
    figure(3)
    if max(stim_mat(:, 4)) < 1
        treatment_trs = find(stim_mat(:, 2) == 1);              %for experiments where some treatment was given to sample (elec stim/changed perfusion etc)
        treatment_odors = unique(stim_mat(treatment_trs, 1));
    else
        treatment_trs = find(stim_mat(:, 4) == 1);
        treatment_odors = unique(stim_mat(treatment_trs, 1));
    end
        
    for odor_ni = 1:length(treatment_odors)
        odor_n = treatment_odors(odor_ni);
        odor_trs = intersect(find(stim_mat(:, 1) == odor_n),  find(stim_mat(:, 2) == 0) );
        odor_treatment_trs = intersect(find(stim_mat(:, 1) == odor_n),  find(stim_mat(:, 2) == 1) );
        
        odor_name = dataset(1).stim.odourNames(odor_n).odour;
        
        subplot(length(treatment_odors), length(treatment_odors), ( (odor_ni-1).*2 + 1))
        imagesc(resp_areas(:, odor_trs), [0, 30])
        
        %averaging response areas across trials within baseline session and
        %listing for all cells
        tr_aved_areas = nanmean(resp_areas(:, odor_trs), 2);
        if isempty(odor_treatment_trs) == 0
            saved_cells_b_shock = [saved_cells_b_shock; tr_aved_areas];
        else
        end
        
        title(['Area under the peak for ' odor_name '; NO-TREATMENT'])
        xlabel('trial number')
        ylabel('cell number')
        try
            subplot(length(treatment_odors), length(treatment_odors), ( (odor_ni-1).*2 + 2) )
        catch
            keyboard
        end
        imagesc(resp_areas(:, odor_treatment_trs), [0, 30])
        
        %averaging response areas across trials within treatment session and
        %listing for all cells
        tr_aved_areas = nanmean(resp_areas(:, odor_treatment_trs), 2);
        if isempty(odor_treatment_trs) == 0
            saved_cells_t_shock = [saved_cells_t_shock; tr_aved_areas];
        else
        end
        
        title(['Area under the peak for ' odor_name '; STIM'])
        xlabel('trial number')
        ylabel('cell number')
        
    end
    set(gcf, 'Color', 'w')
    
    %Plotting averaged response traces for treatment and baseline trials
     figure(5)
    if max(stim_mat(:, 4)) < 1
        treatment_trs = find(stim_mat(:, 2) == 1);              %for experiments where some treatment was given to sample (elec stim/changed perfusion etc)
        treatment_odors = unique(stim_mat(treatment_trs, 1));
    else
        treatment_trs = find(stim_mat(:, 4) == 1);
        treatment_odors = unique(stim_mat(treatment_trs, 1));
    end
        
    for odor_ni = 1:length(treatment_odors)
        odor_n = treatment_odors(odor_ni);
        odor_trs = intersect(find(stim_mat(:, 1) == odor_n),  find(stim_mat(:, 2) == 0) );
        odor_treatment_trs = intersect(find(stim_mat(:, 1) == odor_n),  find(stim_mat(:, 2) == 1) );
        
        odor_name = dataset(1).stim.odourNames(odor_n).odour;
        
        subplot(length(treatment_odors), length(treatment_odors), ( (odor_ni-1).*2 + 1))
        
        baseline_trace_mat = nanmean(dff_data_mat(:, :, odor_trs, odor_n), 3);
        imagesc(baseline_trace_mat', [0, 3]);
        title(['dF/F response traces' odor_name '; NO-TREATMENT'])
        xlabel('trial number')
        ylabel('cell number')
        try
            subplot(length(treatment_odors), length(treatment_odors), ( (odor_ni-1).*2 + 2) )
        catch
            keyboard
        end
        treatment_trace_mat = nanmean(dff_data_mat(:, :, odor_treatment_trs, odor_n), 3);
        imagesc(treatment_trace_mat', [0, 3]);
        title(['dF/F response traces' odor_name '; STIM'])
        xlabel('trial number')
        ylabel('cell number')
        
    end
    set(gcf, 'Color', 'w')
    
    
    
    %plotting significantly responsive cell proportions (sig trace responses on more than half the odor presentations)
    figure(4)
    frac_vec = [];
    for odor_n = 1:n_odors
        sig_cells = sum(sig_cell_mat(:, odor_n))./size(sig_cell_mat, 1);
        
        %skipping if odor wasn't used in expt
        if isnan(sig_cells) == 1
            continue
        else
        end
        
        frac_vec = [frac_vec; sig_cells];
        odor_name = dataset(1).stim.odourNames(odor_n).odour;
        odor_names{1, odor_n} = odor_name;
    end
    unused_odors = isnan(frac_vec);
    
    bar(frac_vec')
    set(gca, 'XtickLabel', odor_names)
    axis([0.5, (length(frac_vec) + 0.5), 0, 1])
    title(['Fraction of sig. responsive cells'])
    xlabel('Odor')
    ylabel('Fraction of responsive cells')

    
    set(gcf, 'Color', 'w')
    try
        switch_trial = dataset(1).stim.switch_valve_trial;
    catch
        switch_trial = [];
    end
    
    if isempty(switch_trial) == 0
        cell_ave_areas = nanmean(resp_areas(:, odor_trs),  2);
        n_cells_b = size(baseline_data_mat, 2);
        
        saved_cells_b = [saved_cells_b; cell_ave_areas(1:n_cells_b)];
        saved_cells_t = [saved_cells_t; cell_ave_areas( (n_cells_b + 1):length(cell_ave_areas) )];
    else
    end
    
    %keyboard
end
fclose(fid);
    
means = [nanmean(saved_cells_b_shock), nanmean(saved_cells_t_shock)];
ses = [nanstd(saved_cells_b_shock)./sqrt(length(saved_cells_b_shock)), nanstd(saved_cells_t_shock)./sqrt(length(saved_cells_t_shock))];

barweb(means', ses')
keyboard


means = [nanmean(saved_cells_b), nanmean(saved_cells_t)];
ses = [nanstd(saved_cells_b)./sqrt(length(saved_cells_b)), nanstd(saved_cells_t)./sqrt(length(saved_cells_t))];

barweb(means', ses')

