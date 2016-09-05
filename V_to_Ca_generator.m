clear all
close all

direc = 'D:\Data\Janelia\Patch\Data_MM\thacq_files\';
list_file = 'cell_list_ABs.xls';
%list_file = 'cell_list_ApBp.xls';
%list_file = 'cell_list_G.xls';
%list_file = 'cell_list_unknown.xls';
%list_file = 'cell_list_unknown_unstained.xls';

odor_list = {'3-Octanol', ...
             '1-Hexanol', ...
             'Pentyl acetate', ...
             '4-Methylcyclohexanol', ...
             '2-Heptanone', ...
             'Diethyl succinate', ...
             'Ethyl lactate', ...
             '1-Octen-3-ol', ...
             'Geranyl acetate', ...
             'Ethyl acetate', ...
             'Empty' ...
            };



[del del1 cell_list] = xlsread([direc list_file]);

n_cells = size(cell_list, 1);
n_odors = size(odor_list, 2);

plot_spike_trace = 0;
plot_single_cell_stuff = 0;

saved_rates = [];
saved_PSTH_curves = [];
saved_sig_cells = [];

for cell_n = 1:n_cells
    cell_path = cell_list{cell_n, 1};
    cell_path = [cell_path, '\'];
    cell_name = cell_list{cell_n, 2};
    
    try

        cell_data = load([cell_path, cell_name, '.mat']);
        spike_data = load([cell_path, cell_name, '_spike.mat']);
        bad_tr_list = load([cell_path, cell_name, '_selection.mat']);
        bad_tr_list = bad_tr_list.rejectedsweeps;
        sf = cell_data.Data_1.parameter.ai_sr;          %AI sampling rate
    catch
        continue
    end
    n_trials = length(fieldnames(cell_data)); 
    
    %initialising stim table, data table
    stim_mat = zeros(n_trials, 3) + nan;        %trial_n x [odor number, odor duration, odor onset time]
    sp_data_mat = zeros(800000, n_trials) + nan;
    sp_wav_mat = [];        
    %loop to load each trial into memory
    for trial_n = 1:n_trials
        
        %skipping trials rejected in THview
        if bad_tr_list(trial_n) == 1
            continue
        else
        end
            
        eval(['curr_data = cell_data.Data_' int2str(trial_n) ';']);
        
        try
            eval(['sp_times = spike_data.Spike_' int2str(trial_n) ';']);
            sp_datapoints = sp_times(:, 1).*sf;
        catch
            sp_times = [];
            sp_datapoints = [];
        end
        
        %adding to stim table
        for odor_ni = 1:n_odors
            list_name = odor_list{1, odor_ni};
            curr_name = curr_data.odor;
            if strcmp(list_name, curr_name) == 1
                curr_odor_n = odor_ni;              %current odor number
                stim_mat(trial_n, 1) = curr_odor_n;
                break
            else
            
            end
        end
        
        stim_mat(trial_n, 2) = curr_data.parameter.odorD;       %odor duration
        stim_mat(trial_n, 3) = curr_data.parameter.preO;        %odor onset time
        stim_mat(trial_n, 4) = curr_data.paratable{2, 3};       %trial duration
        
        v_trace = curr_data.data.voltage;
        sp_vec = zeros((curr_data.parameter.dur.*sf), 1);
        try
            sp_vec(round(sp_datapoints), 1) = 1;
        catch
            keyboard
        end
        sp_data_mat(1:length(sp_vec), trial_n) = sp_vec;
        
        %% Modelling the Ca curve generated by the current V trace
        if curr_data.parameter.odorD < 60
            continue
        else
        end
            
        %CaNew = (1 - fr_time./tau).*Ca + (fr_time./tau)Cbaseline + sp_ht.*n_sp;
        %%from Vogelstein, Paninski et al, J neurophys; DOI: 10.1152/jn.01073.2009
        V_trace = curr_data.data.voltage.*1000/curr_data.amplifier.Vgain;    %in mV
        tau = .5;                   %Ca decay constant in s
        time_step = 1./sf;          %time step in V trace in s
        Cbaseline = .1;             %baseline level of simulated Ca signal
        Ca = Cbaseline;             %initialising Ca signal at t = 0
        n_points = size(V_trace, 1);%n datapoints in current trial
        Vscale = 1./60;             %THIS IS TRICKY; currently mapping -60 mV to 0 mV to 0 to 1 AU. ie 60 mV = 1 AU 
        Ca_max = 1;                 %peak Ca-signal level, dictated by Ca reversal
        
        saved_Ca = zeros(1, n_points) + nan;
        for t_step = 1:n_points
            Vin = V_trace(t_step);
            
            %ensuring Ca-signal is only added by memV, not subtracted by it since VGCCs are not pumps.
            if Vin < -60
                Vin = -60;      %for Vin = -60, scaled added Ca-signal is 0.
            else
            end
                
            CaNew = ( (1 - time_step./tau).*Ca) + ( (time_step./tau).*Cbaseline) + ((Vin + 60).*Vscale);%.*(Ca_max-Ca));
            saved_Ca(1, t_step) = CaNew;
            Ca = CaNew;
        end
        
        figure(1)
        t_vec = curr_data.data.time;
        [hAx,hLine1,hLine2] = plotyy(t_vec, saved_Ca, t_vec, V_trace);

        title('Recorded membrane voltage and simulated Ca-signal')
        xlabel('time (s)')

        ylabel(hAx(1),'simulated Ca-signal (AU)') % left y-axis
        ylabel(hAx(2),'membrane voltage (mV)') % right y-axis
        keyboard
        
    end
end
    n_total_spikes = nansum(nansum(sp_data_mat)); 
    
        
   
    
    