clear all
close all
list_direc = ['C:\Data\CSHL\dataset_list_dopamine_20150421.txt'];
fid = fopen(list_direc);

%loop to go through all experiment datasets listed in list file
while 1
    
    direc = fgetl(fid);
       
    if ischar(direc) ~= 1
        break
    else
    end
    
    direc = [direc '\'];
    %checking if baseline and treatment trials were split into two sets
    if exist([direc 'split_set.txt']) == 2
        split_set = 1;
    else
    end
    
    dataset = load([direc 'expt.mat']);
    dataset = dataset.data;
    raw_data_mat = load([direc 'expt_raw_traces.mat']);
    raw_data_mat = raw_data_mat.raw_data_mat;
    
    
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
        title(['Area under the peak for ' odor_name '; NO-TREATMENT'])
        xlabel('trial number')
        ylabel('cell number')
        try
            subplot(length(treatment_odors), length(treatment_odors), ( (odor_ni-1).*2 + 2) )
        catch
            keyboard
        end
        imagesc(resp_areas(:, odor_treatment_trs), [0, 30])
        title(['Area under the peak for ' odor_name '; STIM'])
        xlabel('trial number')
        ylabel('cell number')
        
    end
    set(gcf, 'Color', 'w')
    
    
    %plotting histograms of significantly responsive cell proportions (sig trace responses on more than half the odor presentations)
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
    
    
    keyboard
end
fclose(fid);
    