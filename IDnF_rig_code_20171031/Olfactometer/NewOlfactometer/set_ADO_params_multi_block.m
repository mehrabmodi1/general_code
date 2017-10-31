function params=set_ADO_params_multi_block()

params=struct;

functionPath=regexp(mfilename('fullpath'),['.*',filesep],'match');
defaults={...
    'duration', [1, 20, 60];...
    'stimLatency',16;...
    'odours',[1, 2, 3, 4, 5];... % for now
    'reps',5;...
    'isi',30;...
    'randomize',1;...
    'firstDilution',0.2;...
    'secondDilution',4500;...
    'post_od_scan_dur', 16;...    
    'n_od_pulses', 1;...                %will deliver n_odor_pulses of length duration and inter_pulse_interval in between
    'inter_pulse_interval', 0.2;...
    'elec_odours', [];...
    'led_odours', [];...
      
    };
    
for P=1:length(params)
    tmpP(P)=addDefaults(params(P));
end
params = tmpP;
for block_n = 2:3
    params = [params, tmpP];
end

%----------------------------------------------------
    function tmp=addDefaults(tmp)
        for i=1:size(defaults,1)
            if ~isfield(tmp,defaults{i,1})
                tmp.(defaults{i,1})=defaults{i,2};
            end
        end              
        %Read in data from the odour list
        [del, tmp.odourNames]=xlsread('E:\Turner lab\Matlab_scripts\Olfactometer\NewOlfactometer\calibration\odorList.xls');
        
    end

end