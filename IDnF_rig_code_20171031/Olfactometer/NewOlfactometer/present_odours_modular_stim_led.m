function present_odours_modular_stim_led(params, scale_isi)

save_dir = curr_aq_direc;
%testing vars
%save_dir = 'D:\Data\testing\olfactometer_code\';

c = onCleanup(@()my_cleanup());        %to shut all valves if user presses Ctrl + C

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
        recovering_session = 1;
    elseif strcmp(button, 'Append') == 1
        
        params_mat_old = load([save_dir 'params.mat']);
        params_mat_old = params_mat_old.params_mat;
       
        %setting up explicit stimulus specification matrix from condensed params structure
        [params_mat_new, params_spec] = setUpStimuli_modular(params);
        params_mat = append_params(params_mat_old, params_mat_new, 0);
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
        [params_mat, params_spec] = setUpStimuli_modular(params);
        start_tr = 1;
    end

else    
    [params_mat, params_spec] = setUpStimuli_modular(params);
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

%reading in manually measured odor propagation delays to align stimuli in time
olf1_olf2_delay = load('E:\Turner lab\Bitbucket_repos\general_code\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\olf1_olf2_delay.mat');
olf1_olf2_delay = olf1_olf2_delay.olf1_olf2_delay;
LED_olf1_delay = load('E:\Turner lab\Bitbucket_repos\general_code\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\LED_olf1_delay.mat');
LED_olf1_delay = LED_olf1_delay.LED_olf1_delay;

for trial_n = start_tr:n_trials
    %reading in stimulus parameters for current trial
    %olfactometer1 params
    odor_n = params_mat(trial_n).odours;
    duration = params_mat(trial_n).duration;
    firstDilution = params_mat(trial_n).firstDilution;
    secondDilution = params_mat(trial_n).secondDilution;
    n_od_pulses = params_mat(trial_n).n_od_pulses;
    inter_pulse_interval = params_mat(trial_n).inter_pulse_interval;
    stim_latency = params_mat(trial_n).stimLatency;
    od_name = params_mat(trial_n).odourNames{odor_n};
    post_od_scan_dur = params_mat(trial_n).post_od_scan_dur;
    pulse_train = params_mat(trial_n).pulse_train;
    pulse_type = params_mat(trial_n).pulse_type;
    stimLatency = params_mat(trial_n).stimLatency;
    
    
    %olfactometer2 params
    if isnan(params_mat(trial_n).duration_olf2) == 1
        no_olf2 = 1;
    elseif isnan(params_mat(trial_n).duration_olf2) == 0
        no_olf2 = 0;
        duration_olf2 = params_mat(trial_n).duration_olf2;
        odor_n_olf2 = params_mat(trial_n).odours_olf2;
        n_od_pulses_olf2 = params_mat(trial_n).n_od_pulses_olf2;
        inter_pulse_interval_olf2 = params_mat(trial_n).inter_pulse_interval_olf2;
        od_name_olf2 = params_mat(trial_n).odourNames_olf2{odor_n_olf2};
        pulse_train_olf2 = params_mat(trial_n).pulse_train_olf2;
        pulse_type_olf2 = params_mat(trial_n).pulse_type_olf2;
        rel_stimLatency_olf2 = params_mat(trial_n).rel_stimLatency_olf2;
    else
    end
    
    
    %listing stimulus parameters and communicating with led/elec stimulating Arduino
    if params_mat(trial_n).led_on == 1
            LED_elec = 0;
            
    elseif params_mat(trial_n).elec_on == 1
        LED_elec = 1;
    else
        LED_elec = 1;
    end
    
    
    if LED_elec == 0
        LED_power = 5.*(params_mat(trial_n).led_power./100);
    elseif LED_elec == 1
        LED_power = 0;
    end
    rel_init_delay = params_mat(trial_n).rel_stim_init_delay;
    init_delay = rel_init_delay.*1000 + stimLatency.*1000 + LED_olf1_delay.*1000;
    duration_ms = params_mat(trial_n).stim_dur.*1000;
    freq_hz = params_mat(trial_n).stim_freq;
    duty_cyc_percent = params_mat(trial_n).st_duty_cyc;

    if isempty(rel_init_delay) == 1
        rel_init_delay = 1000;
        duration_ms = 500;
    else
    end
       
    %communicating stimulus parameters to LED/elec controlling PulsePal
    try
        keyboard
        %FIX PULSEPAL OR REPLACE WITH ARDUINO!!
        %program_pulsepal_LED_elec(LED_elec, init_delay, duration_ms, freq_hz, duty_cyc_percent, LED_power);
    catch
        keyboard
    end
    
    
    %computing total duration of stimulus presentation and imaging for current trial    
    if no_olf2 == 0
        train_dur_olf1 = sum(sum(pulse_train));                                                        %duration of odor train presentation with olf1
        train_dur_olf2 = sum(sum(pulse_train_olf2, 'omitnan'), 'omitnan') + rel_stimLatency_olf2;      %duration of odor train presentation with olf1
        trains_dur = max([train_dur_olf1, train_dur_olf2], [], 'omitnan');                             %duration of longer odor train
        tot_tr_dur = stim_latency + trains_dur + post_od_scan_dur;                                     %total imaging duration (for which scan trigger will be high)
        %computing extra pause needed if olf2 train is longer than olf1 train
        if train_dur_olf2 > train_dur_olf1
            olf2_train_pause = train_dur_olf2 - train_dur_olf1 + post_od_scan_dur;
        else
            olf2_train_pause = 0;
        end
            
    else
        train_dur_olf1 = sum(sum(pulse_train));                                 %duration of odor train presentation with olf1
        trains_dur = train_dur_olf1;
        tot_tr_dur = stim_latency + trains_dur + post_od_scan_dur;          %total imaging duration (for which scan trigger will be high)
        olf2_train_pause = 0;
    end
    
    
    if scale_isi == 0
        isi = params_mat(trial_n).isi;
        if isi < (tot_tr_dur + (od_inj_dur - stimLatency))
            disp('isi is shorter than odor train for olf1 or olf2, make it longer')
            keyboard
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
    
    %olfactometer1 info
    disp(['Trial ' int2str(trial_n) ' of ' int2str(n_trials) '.'])
    disp(['Delivering Odor ' int2str(odor_n) ': ' od_name ' on olfactometer 1.'])
    del_conc = CalcTotalDilution(firstDilution, secondDilution).*100;
    disp(['Concentration delivered ' num2str(del_conc) '%.'])
    if params_mat(trial_n).rand_trains == 0
        disp(['duration ' num2str(duration) 's, n pulses ' int2str(n_od_pulses) '.'])
    elseif params_mat(trial_n).rand_trains == 1
        curr_mean_pulse_dur = params_mat(trial_n).mean_rand_pulse_dur;
        disp(['rand train of train duration ', num2str(duration), ' and mean pulse duration ', num2str(curr_mean_pulse_dur), 's.'])
    else
    end
    
    %olfactometer2 info
    if no_olf2 == 0
        disp(['Delivering Odor ' int2str(odor_n_olf2) ': ' od_name_olf2 ' on olfactometer 2.'])
        if params_mat(trial_n).rand_trains_olf2 == 0
            disp(['duration ' num2str(duration_olf2) 's, n pulses ' int2str(n_od_pulses_olf2) '.'])
        elseif params_mat(trial_n).rand_trains == 1
            curr_mean_pulse_dur = params_mat(trial_n).mean_rand_pulse_dur;
            disp(['rand train of train duration ', num2str(duration_olf2), 's and mean pulse duration ', num2str(curr_mean_pulse_dur), 's.'])
        else
        end
    else
    end
    
       
    
    
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
    addAnalogInputChannel(s,'Dev3', [0, 1], 'Voltage');
    acq_rate = 2000;        %Hz
    s.Rate = acq_rate;
    s.DurationInSeconds = tot_tr_dur;
    lh = addlistener(s,'DataAvailable', @aq_data_bk);
    s.NotifyWhenDataAvailableExceeds = acq_rate.*tot_tr_dur;
    
    
    %communicating stimulus paramters to olfactometer 2 over serial port
    %determining if olf2 needs to be taken out of or put into sleep mode
    if no_olf2 == 0
        if trial_n == 1
            mid_trial = 0;
        elseif trial_n == n_trials
            mid_trial = 2;
        elseif recovering_session == 1
            mid_trial = 0;
        else
            mid_trial = 1;
        end
        odor_vec_olf2 = zeros(size(pulse_train_olf2, 1), 1) + odor_n_olf2;
        initial_delay_olf2 = stimLatency + rel_stimLatency_olf2 - olf1_olf2_delay;
        %sending params to olf2 arduino
        olf_arduino_serial_comm(mid_trial, pulse_train_olf2, odor_vec_olf2, initial_delay_olf2); 
    else
    end
    
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
   
    %pausing for olf2 train to end, if necessary
    pause(olf2_train_pause);
    
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
    recovering_session = 0;
end
release(s)
if exist('PulsePalSystem') == 1
    EndPulsePal;
else
end

%defining clean up function
function [] = my_cleanup()
ShutAllValves_EP;
trigger_scan(0);
if exist('PulsePalSystem') == 1
    EndPulsePal;
else

if no_olf2 == 0
    sleep_olf2              %opens NO valve and closes empty vial valves.
    pause(2);
else
end
    
end
