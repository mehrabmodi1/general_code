function initialiseFlows_MM(AC,firstDilution,secondDilution)
% Set up the MFC flows on the odour delivery system
%
% function initialiseFlows(AC,firstDilution,secondDilution)
%
% Inputs
% AC - the serial port object
% firstDilution - Either a scalar (0 to 1) defining the value for the first 
%                dilution. e.g. 0.1 means 1/10 dilution. Or a vector of
%                length 2, defining flow through the carrier and odour
%                streams. e.g. [0.9,0.3] means 0.9 ml/min carrier and 0.3
%                ml/min odour. 
%                 default.
% second dilution - If NaN, the second dilution is switched off by
%                activating the three 3-port-2-way valves. If not a nan, 
%                this defines the flow rate (in ml/min) which will be added 
%                to the odour flow. By default, this is 5000. 


if nargin<2, firstDilution=0.1; end
if nargin<3, secondDilution=5E3; end

if length(firstDilution)==1
    firstDilution=[1-firstDilution,firstDilution];
end

setFlow(AC,firstDilution(1),'A');
setFlow(AC,firstDilution(2),'B');

%Set flow through carrier
try
    load('E:\Turner lab\Matlab_scripts\Olfactometer\NewOlfactometer\calibration\ADOValvecalib.mat','vCalib');  
    
catch
    fprintf('Can''t load ADOValvecalib.mat, setting carrier to 0.5\n')
    vCalib=0.1; %flow through needle valve 1
    
end

try
    a = load('E:\Turner lab\Matlab_scripts\Olfactometer\NewOlfactometer\calibration\sec_dil_calib.mat');
    secCalib = a.calib_second_dilution;
    
catch
    fprintf('Can''t load sec_dil_calib.mat, setting secondDilution to 4500mLpm\n')
    secCalib=4900; %flow through needle valve 1
end

setFlow(AC,sum(firstDilution)-vCalib,'C')

if secondDilution ~= secCalib
    ans = input('secondDilution being set is different from last calibration. Might need to re-calibrate olfactometer before continuing. Input 0 to abort, anything else to ignore and continue.');
    if ans == 0
        error('aborted for re-calibration!')
    else
    end
else
end

setFlow(AC,secondDilution/5E3,'D') 

% %report dilution
% d1 =firstDilution(2)/sum(firstDilution);
% odourRemaining=(sum(firstDilution)-vCalib)*1000;%this much odour is injected into next stag
% d2=odourRemaining/(secondDilution+odourRemaining);
% totalDilution=d1*d2;
% fprintf('Dilution 1: %0.3g, Dilution 2:%0.3g\nTotal dilution is 1:%4.0f\n',...
%     d1,d2,1/totalDilution)
% 
% if nargout==1    
%     varargout{1}=totalDilution;
% end

    