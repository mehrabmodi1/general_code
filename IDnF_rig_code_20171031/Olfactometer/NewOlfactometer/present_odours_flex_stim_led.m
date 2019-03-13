function present_odours_flex_stim_led(params, scale_isi)

save_dir = curr_aq_direc;
%testing vars
%save_dir = 'D:\Data\testing\olfactometer_code\';

c = onCleanup(@()my_cleanup());        %to shut all valves if user presses Ctrl + C

%recovering interrupted dataset if specified, or starting afresh
%checking if curr direc already has a params file. If yes, prompting user
%to re-specify recover.
if exist([save_dir 'params.mat']) == 2
    button = questdlg('What would you like to do?','Old param file found!','Recover','Append','Over-write', 'Recover');

    if strcmp(button, 'Recover') == 1
        params_mat = load([save_dir 'params.mat']);
        params_mat = params_mat.params_mat;
        n_trials = size(params_mat, 2);
        %identifying last tr completed
        for d_tr_n = 1:n_trials
            done_tr = params_mat(d_tr_n).trs_done;
            if done_tr == 0
                break
            else
            end
        end
        
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
    elseif strcmp(button, 'Append') == 1
        params_mat_old = load([save_dir 'params.mat']);
        params_mat_old = params_mat_old.params_mat;
       
        %setting up explicit stimulus specification matrix from condensed params structure
        [params_mat_new, params_spec] = setUpStimuli_trains_flex(params);
        params_mat = append_params(params_mat_old, params_mat_new);
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
    elseif strcmp(button, 'Over-write') == 1
        [params_mat, params_spec] = setUpStimuli_trains_flex(params);
        start_tr = 1;
    end

else    
    [params_mat, params_spec] = setUpStimuli_trains_flex(params);
    save([save_dir 'params_spec.mat'], 'params_spec');          %saving the params specifications to file
    start_tr = 1;
end

n_trials = size(params_mat, 2);

ShutAllValves_EP;       %making sure all odor valves are closed and carrier stream is directed through shuttle valve to fly.

%opening serial port to communicate with the Alicats
AC = instrfind('Tag','Alicat_serial');  %checking if Alicat Serial port is already open.
if isempty(AC) == 1
    AC = connectAlicat();
end
    
secondDilution1 = params_mat(1).secondDilution;
initialiseFlows_MM(AC, 0.1, secondDilution1);      %initialising flows for the first time just to set things up. 
od_inj_dur = 24.5;                         %this is the duration in s for which MFC B flow is injected into an odor vial to fully fill the system with odor. Stim_latency has to be longer than this.

for trial_n = start_tr:n_trials
  
    odor_n = params_mat(trial_n).odours;
    duration = params_mat(trial_n).duration;
    firstDilution = params_mat(trial_n).firstDilution;
    secondDilution = params_mat(trial_n).secondDilution;
    n_od_pulses = params_mat(trial_n).n_od_pulses;
    inter_pulse_interval = params_mat(trial_n).inter_pulse_interval;
    stim_latency = params_mat(trial_n).stimLatency;    
    od_name = params_mat(trial_n).odourNames{odor_n};
    post_od_scan_dur = params_mat(trial_n).post_od_scan_dur;
    pulse_train = params_mat(trial_n).rand_train;
    pulse_type = params_mat(trial_n).pulse_type;
    LED_power = 5.*(params_mat(trial_n).led_power./100);
    
    
    %listing parameters and communicating with led/elec stimulating Arduino
    %if necessary for current trial
    if params_mat(trial_n).led_on == 1 || params_mat(trial_n).elec_on == 1      %condition to check if LED or elec stim is being delivered on current trial
        if params_mat(trial_n).led_on == 1
            LED_elec = 0;
            
        elseif params_mat(trial_n).elec_on == 1
            LED_elec = 1;
        end
        init_delay_ms = params_mat(trial_n).stim_init_delay_ms;
        duration_ms = params_mat(trial_n).stim_dur;
        freq_hz = params_mat(trial_n).stim_freq;
        duty_cyc_percent = params_mat(trial_n).st_duty_cyc;
        
        %communicating stimulus parameters to Arduino
        keyboard
        stim_arduino_serial_comm(LED_elec, init_delay_ms, duration_ms, freq_hz, duty_cyc_percent, LED_power);
        
            
    else
    end
    
    if pulse_type == 1 
        tot_tr_dur = stim_latency + duration + ...
            (duration + inter_pulse_interval).*(n_od_pulses - 1) + post_od_scan_dur;
    elseif pulse_type == 0
        tot_tr_dur = stim_latency + duration + ...
            (duration + duration).*(n_od_pulses - 1) + post_od_scan_dur;
    else
    end
    
    if scale_isi == 0
        isi = params_mat(trial_n).isi;
    elseif scale_isi == 1
        isi = max([60, ((tot_tr_dur - stim_latency - post_od_scan_dur).*3)]);         %scales isi to stim duration, with a minimum isi of 60s
        isi = min([isi, 130]);                                                        %making sure scaled isi doesn't exceed 130s
        params_mat(trial_n).isi = isi;
    end
    
    %displaying total time until end of acqn
    n_trials_left = n_trials - trial_n + 1;
    tot_time = round((n_trials_left.*isi)./60);  %in min
    disp([int2str(tot_time), ' minutes left for completion of acquisition.']);
    
    
    disp(['Trial ' int2str(trial_n) ' of ' int2str(n_trials) '.'])
    disp(['Delivering Odor ' int2str(odor_n) ': ' od_name '.'])
    del_conc = CalcTotalDilution(firstDilution, secondDilution).*100;
    disp(['Concentration delivered ' num2str(del_conc) '%.'])
    disp(['duration ' num2str(duration) 's, n pulses ' int2str(n_od_pulses) '.'])
    
    
    %checking to see if olfactometer odor fill needs to be initialised
    %before scan trigger or after.
    if stim_latency < od_inj_dur
        od_fill_early = 1;
    elseif stim_latency >= od_inj_dur
        od_fill_early = 0;
    else
    end


    %% delivering odor
    %Setting up PID acuisition, 
    s = daq.createSession('ni');
    addAnalogInputChannel(s,'Dev3', 0, 'Voltage');
    acq_rate = 2000;        %Hz
    s.Rate = acq_rate;
    s.DurationInSeconds = tot_tr_dur;
    lh = addlistener(s,'DataAvailable', @aq_data_bk);
    s.NotifyWhenDataAvailableExceeds = acq_rate.*tot_tr_dur;
    
    initialiseFlows_MM(AC, firstDilution, secondDilution);  %setting MFC flow rates for required conc.
    
    if od_fill_early == 0
        tic
        t_stamp = now;
        s.startBackground();                    %starting PID acqn in the background
        trigger_scan(1);                        %triggering ScanImage to start image acquisition
        disp('pre odor scanning...')
        pause(stim_latency - od_inj_dur)        %pause before filling system with odor for long stim latencies
        injectOdour_EP(odor_n)                  %filling system with odor, switching MFC B flow from empty vial to odor vial

        pause(od_inj_dur)                       %waiting for system to get filled with odor
    
    elseif od_fill_early == 1
        injectOdour_EP(odor_n)                  %filling system with odor, switching MFC B flow from empty vial to odor vial
        pause((od_inj_dur - stim_latency));     %waiting for olfactometer odorisation before triggering a shorter stim_latency scan
        tic
        t_stamp = now;
        s.startBackground();                    %starting PID acqn in the background
        trigger_scan(1);                        %triggering ScanImage to start image acquisition
        disp('pre odor scanning...')
        pause(stim_latency)   
    else
    end
    
    %flipping shuttle valve to deliver odor pulse(s)
    disp('odor being delivered...')
    for r_pulse_n = 1:size(pulse_train, 1)
        pause(pulse_train(r_pulse_n, 1));
        FlipValve_EP('Final',0)        
        
        pause(pulse_train(r_pulse_n, 2));
        FlipValve_EP('Final',1)
        
    end
    
       
    ShutAllValves_EP;
    disp('post odor scanning...')
    pause(post_od_scan_dur)                 %waiting to end image acquisition
    trigger_scan(0)                         %ending image acquisition
    
       
    %logging current trial as done and saving params_mat
    params_mat(trial_n).trs_done = t_stamp;     %time stamp recorded at the beginning of the trial
    save([save_dir 'params.mat'], 'params_mat');                %saving the detailed parameters for each trial to file
    
    
    
    %re-naming PID trace file saved in the background
    PID_data = load([save_dir 'last_PID_trace.mat']);
    PID_data = PID_data.data;
    save([save_dir 'PID_trace_tr-' int2str(trial_n) '.mat'], 'PID_data');
    
    disp('Updated param-file. Waiting for isi.');
    disp(' ')
    disp(' ')
    
    %tricking scanimage into releasing current trial file...
    %by triggering a fake, short trial
    pause(1)
    trigger_scan(1)
    pause(1)
    trigger_scan(0)    
    
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
    
end
release(s)

%defining clean up function
function [] = my_cleanup()
ShutAllValves_EP;
trigger_scan(0);
