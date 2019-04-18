load ParameterMatrix_Example; % Loads the default parameter matrix
ProgramPulsePal(ParameterMatrix); % Sends the default parameter matrix to Pulse Pal
ProgramPulsePalParam(1, 'Phase1Voltage', 2.5); % Set output channel 1 to produce 2.5V pulses
ProgramPulsePalParam(1, 'Phase1Duration', 0.002); % Set output channel 1 to produce 2ms pulses
ProgramPulsePalParam(1, 'InterPulseInterval', 0.098); % Set pulse interval to produce 10Hz pulses
ProgramPulsePalParam(1, 'PulseTrainDuration', 120); % Set pulse train to last 120 seconds
ProgramPulsePalParam(1, 12, 1); % Set output channel 1 to respond to trigger ch 1
ProgramPulsePalParam(1, 'TriggerMode', 1); % Set trigger channel 1 to toggle mode