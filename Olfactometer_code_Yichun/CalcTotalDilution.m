function TotalDilution=CalcTotalDilution(firstdilution,seconddilution)
% Report total dilution.
%
% firstDilution - Either a scalar (0 to 1) defining the value for the first 
%                dilution. e.g. 0.1 means 1/10 dilution. Or a vector of
%                length 2, defining flow through the carrier and odour
%                streams. e.g. [0.9,0.3] means 0.9 ml/min carrier and 0.3
%                ml/min odour. 
%
% second dilution - Defines the flow rate (in ml/min) which will be added 
%                to the odour flow.
% 

try
    load('C:\Users\Glenn\Documents\MATLAB\ADOValvecalib.mat','vCalib');    
catch
    fprintf('Can''t load ADOValvecalib.mat, setting carrier to 0.05\n')
    vCalib=0.95; %flow through needle valve 1
end

if length(firstdilution)==1
    firstdilution=[1-firstdilution,firstdilution];
end

d1 =firstdilution(2)/sum(firstdilution);
odourRemaining=(sum(firstdilution)-vCalib)*1000;%this much odour is injected into next stag
d2=odourRemaining/(seconddilution+odourRemaining);
TotalDilution=d1*d2;