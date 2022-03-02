function stripped_down_present_stim(params, scale_isi)

save_dir = curr_aq_direc;
%testing vars
%save_dir = 'D:\Data\testing\olfactometer_code\';

c = onCleanup(@()my_cleanup());        %Things to ahut down and clean up when Ctrl + C

%recovering interrupted dataset if specified, or starting afresh
%checking if curr direc already has a params file. If yes, prompting user
%to re-specify recover.
recovering_session = 0;
if exist([save_dir 'params.mat']) == 2
    button = questdlg('What would you like to do?','Old param file found!','Recover','Append','Over-write', 'Recover');

    if strcmp(button, 'Recover') == 1
        params_mat = load([save_dir 'params.mat']);
        params_mat = params_mat.params_mat;
        n_trials = size(params_mat, 2);
        %identifying last tr completed
        for d_tr_n = n_trials:-1:1
            done_tr = params_mat(d_tr_n).trs_done;
            if done_tr ~= 0
                broken_loop = 1;
                break
            else
            end
            broken_loop = 0;
        end 
        
        if broken_loop == 1
            d_tr_n = d_tr_n + 1;
        else
        end
        
        %asking user to specify trial n to recover session from
        a = inputdlg('Input trial n to acquire next; 0 for last trial done.', 'Trial n', 1, {'0'});
        a = str2num(a{1, 1});
        if a == 0
            start_tr = d_tr_n;
        else
            start_tr = a;
        end
        
        for tr_n = start_tr:n_trials
            params_mat(tr_n).trs_done = 0;
        end
        recovering_session = 1;
    elseif strcmp(button, 'Append') == 1        %appending trials specified by params variable passed to the master program function to the end of the recovered trial set
        
        params_mat_old = load([save_dir 'params.mat']);
        params_mat_old = params_mat_old.params_mat;
       
        [params_mat_new, params_spec] = setUpStimuli_modular(params);       %setting up explicit, trial by trial stimulus specifications for new params variable passed to master function
        params_mat = append_params(params_mat_old, params_mat_new, 0);      %write short function to append parameter structures depending on how they are specified
        n_trials = size(params_mat, 2);
        
        %identifying last tr completed
        for d_tr_n = 1:n_trials
            done_tr = params_mat(d_tr_n).trs_done;
            if done_tr == 0
                break
            else
            end
        end
        start_tr = d_tr_n;
        save([save_dir 'params_spec2.mat'], 'params_spec');          %saving the params specifications to file
        recovering_session = 1;
    elseif strcmp(button, 'Over-write') == 1
        [params_mat, params_spec] = setUpStimuli_modular(params);   %setting up explicit, trial by trial stimulus specifications
        start_tr = 1;
    end

else    
    [params_mat, params_spec] = setUpStimuli_modular(params);
    save([save_dir 'params_spec.mat'], 'params_spec');          %saving the params specifications to file
    start_tr = 1;
end

n_trials = size(params_mat, 2);


for trial_n = start_tr:n_trials
    %reading in stimulus parameters for current trial
    post_od_scan_dur = params_mat(trial_n).post_od_scan_dur;
    stimLatency = params_mat(trial_n).stimLatency;
    trig_scan = params_mat(trial_n).trigger_scan;
    
    %listing stimulus parameters and communicating with led/elec stimulating Arduino
    if params_mat(trial_n).led_on == 1
            LED_elec = 0;
            
    elseif params_mat(trial_n).elec_on == 1
        LED_elec = 1;
    else
        LED_elec = 2; %niether LED nor elec stim to be delivered
    end
    
    
    if LED_elec == 0
        LED_power = 5.*(params_mat(trial_n).led_power./100);
    elseif LED_elec == 1
        LED_power = 0;
    end
    rel_init_delay = params_mat(trial_n).rel_stim_init_delay;
    init_delay_ms = rel_init_delay.*1000 + stimLatency.*1000 + LED_olf1_delay.*1000;        %This will need to be re-specified
    duration_ms = params_mat(trial_n).stim_dur.*1000;    
    freq_hz = params_mat(trial_n).stim_freq;
    duty_cyc_percent = params_mat(trial_n).st_duty_cyc;

    if isempty(rel_init_delay) == 1
        rel_init_delay = 1000;
        duration_ms = 500;
    else
    end
       
    %communicating stimulus parameters to LED/elec controlling arduino
    disp('warning: The LED_power param doesn''t actually control LED power. This is currently adjusted manually to 5% with a V-divider.')
    stim_arduino_serial_comm(LED_elec, init_delay_ms, duration_ms, freq_hz, duty_cyc_percent);
     
    
    
    %computing total duration of stimulus presentation and imaging for current trial    
    %user needs to write code here to compute total trial duration
    tot_time = user_parameter1 + user_parameter2 + ;   %in minutes       
    
    %code to automatically re-scale inter trial interval depending on the total stimulus duration on a given trial
    if scale_isi == 0
        isi = params_mat(trial_n).isi;
        if isi < (tot_tr_dur + (od_inj_dur - stimLatency))
            disp('isi is shorter than odor train for olf1 or olf2, make it longer')
            isi = (tot_tr_dur + (od_inj_dur - stimLatency)) + 2;
            params_mat(trial_n).isi = isi;
            
        else
        end
        
    elseif scale_isi == 1
        isi = max([60, ((tot_tr_dur - stim_latency - post_od_scan_dur).*3)]);         %scales isi to stim duration, with a minimum isi of 60s
        isi = min([isi, 130]);                                                        %making sure scaled isi doesn't exceed 130s
        params_mat(trial_n).isi = isi;
    end
    
    
    
    %displaying total time until end of acqn and other trial info
    n_trials_left = n_trials - trial_n + 1;
    tot_time = round((n_trials_left.*isi)./60);  %in min
    disp(['approx.', int2str(tot_time), ' minutes left for completion of acquisition.']);
    
        
    

    %% delivering stimuli
    %Setting up analog signal acuisition, 
    s = daq.createSession('ni');
    ch = addAnalogInputChannel(s,'Dev3', [0, 2, 3], 'Voltage');        %User will have to specify source of analog signal (NI DAQ used here)
    %ch = addAnalogInputChannel(s,'PXI2Slot2_chs2', [2, 3], 'Voltage');
    ch(1).Coupling = 'DC';
    
    acq_rate = 2000;        %Hz
    s.Rate = acq_rate;
    s.DurationInSeconds = tot_tr_dur;
    lh = addlistener(s,'DataAvailable', @aq_data_bk);
    s.NotifyWhenDataAvailableExceeds = acq_rate.*tot_tr_dur;
    
    s.startBackground();                    %starting PID acqn in the background
        
    if trig_scan == 1
        trigger_scan(1, 1);                        %triggering ScanImage to start image acquisition
        disp('pre odor scanning...')
    else
        trigger_scan(0, 1);
    end
    
    %Insert code here to execute anything else during trial period
    
    if trig_scan == 1
        disp('post odor scanning...')
        pause(post_od_scan_dur)                 %waiting to end image acquisition
    else
        pause(10)  
    end
    trigger_scan(0, 0);         %ending image acquisition
    
    
       
    %logging current trial as done and re-writing params_mat to file
    params_mat(trial_n).trs_done = t_stamp;     %time stamp recorded at the beginning of the trial
    save([save_dir 'params.mat'], 'params_mat');                %saving the detailed parameters for each trial to file
    
    
    %re-naming analog trace file saved in the background
    PID_data = load([save_dir 'last_PID_trace.mat']);
    PID_data = PID_data.data;
    save([save_dir 'PID_trace_tr-' int2str(trial_n) '.mat'], 'PID_data');
    
    disp('Updated param-file. Waiting for isi.');
    disp(' ')
    disp(' ')
    
    %tricking scanimage into releasing current trial file...
    %by triggering a fake, short trial
    if trig_scan == 1
        pause(1)
        trigger_scan(1, 0)
        pause(1)
        trigger_scan(0, 0)    
    else
        pause(2)
    end
    
    
    %pause for inter stimulus interval (between this and next trial)
    if trial_n < n_trials
        if od_fill_early == 0
            pause(isi-toc)
        elseif od_fill_early == 1
            pause(isi - toc - (od_inj_dur - stim_latency))      %accounting for the extra seconds spent odorising the olfactometer before triggering scan.
        else
        end
    else
    end
    recovering_session = 0;
end
pause(4)
release(s)
close_serial_port(19);   %LED_arduino


%defining clean up function
function [] = my_cleanup()

trigger_scan(0, 0);
close_serial_port(19)   %LED_arduino
    
