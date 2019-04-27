function [] = program_pulsepal_LED_elec(LED_elec, init_delay_ms, duration_ms, freq_hz, duty_cyc_percent, stim_voltage)

%connecting to pulse pal if not already connected
if exist('PulsePalSystem') == 0
    PulsePal;
else
end

load ParameterMatrix_Example; % Loads the default parameter matrix
ProgramPulsePal(ParameterMatrix); % Sends the default parameter matrix to Pulse Pal
LED_warning_lead_t = 0.02; %(in s)

%specifying stim_train channel number
if LED_elec == 0
    channel_on = 1;
    channel_off = 2;
elseif LED_elec == 1
    channel_on = 2;
    channel_off = 1;
else
end

%computing pulse train parameters
cyc_dur = 1./freq_hz;
on_dur = cyc_dur.*(duty_cyc_percent./100);
off_dur = cyc_dur - on_dur;
init_delay = init_delay_ms./1000;
duration = duration_ms./1000;

%computing LED warning signal parameters to control PMT shutters
if off_dur > LED_warning_lead_t
    on_dur_w = on_dur + (2.*LED_warning_lead_t);
    off_dur_w = off_dur - (2.*LED_warning_lead_t);
    init_delay_w = init_delay - LED_warning_lead_t;
elseif off_dur <= LED_warning_lead_t                %case where there isn't enough inter-pulse-interval to generate a shutter warning
    on_dur_w = duration + 2.*LED_warning_lead_t;
    off_dur_w = 0;
    init_delay_w = init_delay - LED_warning_lead_t;
else
end

%only triggering LED_shutter warning signal if an LED stimulus is being
%delivered.
if LED_elec == 0
    warning_V = 5;
else
    warning_V = 0;
end
    
%rounding to nearest 100 microseconds
on_dur = round(on_dur, 4);
off_dur = round(off_dur, 4);
duration = round(duration, 4);
init_delay = round(init_delay, 4);
on_dur_w = round(on_dur_w, 4);
off_dur_w = round(off_dur_w, 4);
init_delay_w = round(init_delay_w, 4);

%communicating pulse train parameters to PulsePal for stimulus on channel
bitout = ProgramPulsePalParam(channel_on, 'Phase1Voltage', stim_voltage); % Set
bitout = ProgramPulsePalParam(channel_on, 'Phase1Duration', on_dur); % Set output channel 1 to produce pulses
bitout = ProgramPulsePalParam(channel_on, 'InterPulseInterval', off_dur); % Set inter pulse interval to off_dur
bitout = ProgramPulsePalParam(channel_on, 'PulseTrainDuration', duration); % Set pulse train to last duration seconds
bitout = ProgramPulsePalParam(channel_on, 'PulseTrainDelay', init_delay); % Set pulse train onset delay
bitout = ProgramPulsePalParam(channel_on, 12, 1); % Set output channel to respond to trigger ch 1
bitout = ProgramPulsePalParam(1, 'TriggerMode', 2); % Set trigger channel to pulse-gated mode

%communicating pulse train parameters to PulsePal for stimulus off channel
bitout = ProgramPulsePalParam(channel_off, 'Phase1Voltage', 0); % Set
bitout = ProgramPulsePalParam(channel_off, 'Phase1Duration', on_dur); % Set output channel 1 to produce 2ms pulses
bitout = ProgramPulsePalParam(channel_off, 'InterPulseInterval', off_dur); % Set inter pulse interval to off_dur
bitout = ProgramPulsePalParam(channel_off, 'PulseTrainDuration', duration); % Set pulse train to last duration seconds
bitout = ProgramPulsePalParam(channel_off, 'PulseTrainDelay', init_delay); % Set pulse train onset delay
bitout = ProgramPulsePalParam(channel_off, 12, 1); % Set output channel to respond to trigger ch 1


%communicating pulse train parameters to PulsePal for LED warning signal to shutter control arduino
bitout = ProgramPulsePalParam(3, 'Phase1Voltage', warning_V); % Set
bitout = ProgramPulsePalParam(3, 'Phase1Duration', on_dur_w); % Set output channel 1 to produce 2ms pulses
bitout = ProgramPulsePalParam(3, 'InterPulseInterval', off_dur_w); % Set inter pulse interval to off_dur
bitout = ProgramPulsePalParam(3, 'PulseTrainDuration', duration); % Set pulse train to last duration seconds
bitout = ProgramPulsePalParam(3, 'PulseTrainDelay', init_delay_w); % Set pulse train onset delay
bitout = ProgramPulsePalParam(3, 12, 1); % Set output channel to respond to trigger ch 1

