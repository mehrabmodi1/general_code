function params=set_ADO_params_find_resp_ods(params)
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
% Rob Campbll - July 2010

  
%Set up default parameters
if nargin==0
  params=struct;
end


functionPath=regexp(mfilename('fullpath'),['.*',filesep],'match');
defaults={...
    'duration', [0.2, 1, 5, 10];...
    'stimLatency',25;...
    'odours',[3, 5];... % for now
    'reps',1;...
    'isi',55;...
    'randomize',1;...
    'firstDilution',0.2;...
    'secondDilution',4500;...
    'post_od_scan_dur', 7;...    
    'n_od_pulses', 5;...                %will deliver n_odor_pulses of length duration and inter_pulse_interval in between
    'inter_pulse_interval', 1.5;...
    'elec_odours', [];...
    'led_odours', [];...
    'rand_trains', 0;
    'n_rand_trains', 1;
    'min_pulse_dur', 0.1;           %in seconds, the shortest pulse the olfactometer can deliver
    'pulse_type', 0;                %0 - duty cycle of 0.5, ignores inter_pulse_interval; 1 - use inter_pulse_interval
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
        
    end

end