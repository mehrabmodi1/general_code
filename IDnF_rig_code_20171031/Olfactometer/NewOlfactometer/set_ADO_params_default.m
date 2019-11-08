function params=set_ADO_params_default(params)
%  function params=defaultADOparams
%
% Create a parameter structure with the default values specified by
% the body of this function. Any of these parameters can be
% modified or others added providing a partial parameters structure
% as an argument to this function. 
%
% e.g.
% params.odours=[6,9,14];
% params.reps=10;  
% params=defaultADOparams(params)
%
% The params object is then fed to ado_deliverOdours, to define what will 
% be presented to the fly. Note that within that function we convert the
% structure to stimulus matrix using this:
% params=setUpStimuli(params);
%
% Rob Campbell - July 2010, re-written by Mehrab Modi - 2016, 2018, 2019

  
%Set up default parameters
if nargin==0
  params=struct;
end


functionPath=regexp(mfilename('fullpath'),['.*',filesep],'match');
defaults={...
    
    %common params
    'reps',6;...
    'isi',50;...
    'randomize',1;...
    'post_od_scan_dur', 10;... 
    'stimLatency', 10;...
    
    %olfactometer1 params
    'duration', [0.05];...                %for rand trains - duration of train, for simple trains, duration of each pulse
    'odours',[3, 4];...                    %olf1_list = [3, 5, 10, 11] corress olf2_list = [1, 3, 2, 4]
    'firstDilution',0.1;...
    'secondDilution',4600;...
    %simple train params
    'n_od_pulses', 1;...                %will deliver n_odor_pulses of length duration and inter_pulse_interval in between
    'inter_pulse_interval', 1.5;...
    'pulse_type', 0;                    %0 - duty cycle of 0.5, ignores inter_pulse_interval; 1 - use inter_pulse_interval
    %random train params
    'rand_trains', 0;
    'n_rand_trains', 2;
    'min_pulse_dur', 0.1;               %in seconds, the shortest pulse the olfactometer can deliver
    'mean_rand_pulse_dur', 2;                %mean pulse duration in rand trains for olfactometer 1
    
        
    %olfactometer2 params
    'duration_olf2', [3];...            %for rand trains - duration of train, for simple trains, duration of each pulse
    'odours_olf2',[2, 3];...                %olf1_list = [3, 5, 10, 11] corress olf2_list = [1, 3, 2, 4]
    'rel_stimLatency_olf2', 2;...       %stim latency of olfactometer2 pulse train relative to olfactometer1 pulse train
    %simple train params
    'n_od_pulses_olf2', 1;...           %will deliver n_odor_pulses of length duration and inter_pulse_interval in between
    'inter_pulse_interval_olf2', 5;...
    'pulse_type_olf2', 1;               %0 - duty cycle of 0.5, ignores inter_pulse_interval; 1 - use inter_pulse_interval
    %rand_train params
    'rand_trains_olf2', 0;
    'n_rand_trains_olf2', 2;            %number of different random train instantiations
    'min_pulse_dur_olf2', 0.5;          %in seconds, the shortest pulse the olfactometer can deliver
    'mean_rand_pulse_dur_olf2', 2;      %mean pulse duration in rand trains for olfactometer 2
    
    %led or elec stim params
    'elec_odours', [];...
    'led_odours', [];...
    'LED_power', [5];
    'rel_stim_init_delay', [0];...
    'stim_dur', [10];...
    'stim_freq', 1;...
    'st_duty_cyc', 50;...
    
    
    };
    


for P=1:length(params)
    tmpP(P)=addDefaults(params(P));
end
params=tmpP;

%----------------------------------------------------
    function tmp=addDefaults(tmp)
        for i=1:size(defaults,1)
            if ~isfield(tmp,defaults{i,1})
                tmp.(defaults{i,1})=defaults{i,2};
            end
        end
        
        %Read in data from the odour list
        fid = fopen('odor_list_path.txt');
        odor_list_path = fgetl(fid);
        fclose(fid);
        [del, tmp.odourNames]=xlsread(odor_list_path);
        
        %Read in data from the odour list for olf2
        fid = fopen('odor_list_path_olf2.txt');
        odor_list_path_olf2 = fgetl(fid);
        fclose(fid);
        [del, tmp.odourNames_olf2]=xlsread(odor_list_path_olf2);
        
    end

end