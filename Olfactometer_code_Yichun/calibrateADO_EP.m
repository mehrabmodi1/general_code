function calibrateADO_EP(AC)
% Modified from calibrateADO.m 
% Alternate between settings required for odour delivery or for calibrating the first needle vavle
%
% function calibrateNeedleValve1(AC)
%
% Inputs
% AC - the serial port object


%Read information from MFCs to store initial states
fprintf(AC,'A'); A=readMFC(AC);
fprintf(AC,'B'); B=readMFC(AC);
setFlow(AC,0,'AB')

%Now we set MFC A to close loop pressure
setPressure(AC,'A',14.71)
%and stop air from escaping through the air vial and exhast
FlipValve_EP('NO',0) %close the NO clean air
%block exhaust so all flow is through needle valve
FlipValve_EP('Block',0)

%The user now adjusts the needle valve
fprintf(['Adjust second dilution needle valve\n',...
    'whilst reading volumetric flow on MFC A.\n',...
    'Then press RETURN.\n'])
pause

%when it's oscilating a lot we can't tak samples and it's not stable. So
%only start to record once 5 seconds of good data have been obtained. 
good=0;
step=0.25;
checkLim=5;
fprintf('Checking that reading is stable\n')
while good<checkLim
    fprintf(AC,'A');
    tmp=readMFC(AC);
    if isempty(tmp),
        pause(step)
        good=0;
        checkLim=8; %Become cautious if a problem is found
        continue,end
    good=good+step;
    pause(step)
end


%We now measure the mean volumetric flow through this valve
sampleDuration=30; %seconds
samplesPerSecond=2;
tic
data=nan(1,sampleDuration*samplesPerSecond);
n=1;
clf
X=@(L) (0:1/samplesPerSecond:(L/samplesPerSecond)-1/samplesPerSecond);
p=plot(X(length(data)),data,'-k.');
xlim([0,sampleDuration])
ylim([0.2,1.1])

while toc<=sampleDuration+1/samplesPerSecond
    fprintf(AC,'A');
    tmp=readMFC(AC);
    if isempty(tmp),continue,end
    data(n)=tmp.volumetricFlow;
    n=n+1;
    set(p,'YData',data,'XData',X(length(data)))
    pause(1/samplesPerSecond);
end


ylim([min(data)-min(data)*0.25,...
    max(data)+max(data)*0.15])
    
hold on
mu=nanmean(data);
plot(xlim,[mu,mu],'r--','linewidth',2)
hold off

fprintf('Calibrated! Stored mean needle valve flow of %0.3f ml/min\n',mu)
vCalib=mu;
save('C:\Users\Glenn\Documents\MATLAB\thacq.ver4\ADOValvecalib.mat','vCalib')


%Open valves and return flows to where the were
FlipValve_EP('NO',1) %close the NO clean air
%block exhaust so all flow is through needle valve
FlipValve_EP('Block',1)
writeLoop(AC,'A',2) %back to volumetric flow
setFlow(AC,A.volumetricFlow,'A')
setFlow(AC,B.volumetricFlow,'B')