function params=defaultADOparamsMehrab_stim_dur_20160303(params)
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
params=tmpP;

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