function [temp,fitted] = IsItLinear(cellinfo,mode,plotmode)
if nargin<3; plotmode = 'off'; end
% defaults
binsize = 50; % in msec
mode = mode - 1; % odor A, B or C
protoIDs = [(2+mode) (9+mode) (21+mode) 24 25 26]; 
count = 0;

% getting all the needed data from the spike files and odor time stamps
% enter the main data directory
cd (cellinfo.path);
% enter the data directory for the cell in concern
cd ([cellinfo.recdate,'/processed/',cellinfo.datafile]);

if any(size(dir([pwd '/*.params' ]),1))
    [dataALL,subtypesALL,odorIDs,list,settings] = ...
    MTanalyze2014(cellinfo.tetrode,cellinfo.ID,'all',plotmode);  
else
    [dataALL,subtypesALL,odorIDs,list,settings] = ...
    MTanalyze2013(cellinfo.tetrode,cellinfo.ID,'all',plotmode); 
end
Wvfrmmode = 1;

startveclength = 1000*(0.05 + settings.suction(1))/binsize;
protocols = fieldnames(dataALL);

% for each protocol, get the structure temp that has all 
% required data and stimulus waveforms
for proto = 1:length(protoIDs)
    % find the protocol name in dataALL
    protoname = (char(protocols(find(list(:,2)==protoIDs(proto)))));
if ~isempty(protoname)
    if proto>=4
        protonickname = [protoname(1),protoname(4)];
    else
        protonickname = protoname(1:2);
    end
        
    data = dataALL.(protoname);  
    % this gives me all the subtypes of the protocol
    % for each subtype
    for subtype = 1:length(fieldnames(data))
        if ~((proto>=4)&&(isemptyrandom(mode,subtypesALL.(protoname)(subtype,1))))
        raster = data.(['f',num2str(subtype)]).raster;
        plotwindow = data.(['f',num2str(subtype)]).plotwindow;
        % convert each raster
        [FR] = raster2FR(raster,plotwindow,binsize);
        %FR.raw.mean,FR.raw.sd,FR.raw.sem,FR.raw.noise,FR.raw.bins
        %FR.smooth.mean,FR.smooth.sd,FR.smooth.sem,FR.smooth.noise,
        %FR.smooth.bins
        bins = FR.raw.bins;        
        
        count = count + 1;
        temp.raw.mean(count,1:bins) = FR.raw.mean;
        temp.raw.sd(count,1:bins) = FR.raw.sd;
        temp.raw.sem(count,1:bins) = FR.raw.sem;
        temp.raw.noise(count,1:bins) = FR.raw.noise;
        
        bins = FR.smooth.bins;
        temp.smooth.mean(count,1:bins) = FR.smooth.mean;
        temp.smooth.sd(count,1:bins) = FR.smooth.sd;
        temp.smooth.sem(count,1:bins) = FR.smooth.sem;
        temp.smooth.noise(count,1:bins) = FR.smooth.noise;
        
        % which protocol and which variable
        temp.protocol(count,:) = protonickname;
        temp.protocolID(count,1) = protoIDs(proto);
        temp.variable(count,1) = subtypesALL.(protoname)(subtype,1);
        
        % odorwaveform 
        [wfrm] = OdorWaveform(protonickname,subtypesALL.(protoname)(subtype,1),...
        settings.suction,mode,...
        Wvfrmmode,binsize,odorIDs(mode+1));
        temp.odorwfm(count,1:length(wfrm)) = wfrm;
        temp.odorlength(count,1) = length(wfrm);
        
        % start vector
        if (proto == 1)&&(subtypesALL.(protoname)(subtype,1)==0.2) % single pulse, 200ms
            temp.raw.startvector = FR.raw.mean(1,1:startveclength);
            temp.smooth.startvector = FR.smooth.mean(1,1:startveclength);
        end
        end
    end
end
end

temp.baseline = settings.baseline;

%% use error minimization to find the kernel
fitmode = 'smooth';

% assemble startvector if not already allocated
if ~isfield(temp.(fitmode),'startvector')
    temp.(fitmode).startvector = zeros(1,startveclength);
end

% assemble xdata, ydata for lsqcurvefit
fprintf (1, '\nfitting ....\n');
[fitted] = ForGetKernel(temp,fitmode);

% use all protocols
[fitted.kernel1.kernel] = GetKernel(temp.(fitmode).startvector,temp.baseline,fitted.xdata,fitted.ydata);

% use only random and pulse trains
f = find(temp.protocolID>20);
if ~isempty(f)
[fitted.kernel2.kernel] = GetKernel(temp.(fitmode).startvector,temp.baseline,fitted.xdata(:,f),fitted.ydata(:,f));
end

% all but random and pulse trains
f = find(temp.protocolID<20);
if ~isempty(f)
[fitted.kernel3.kernel] = GetKernel(temp.(fitmode).startvector,temp.baseline,fitted.xdata(:,f),fitted.ydata(:,f));
end

% fitted.zdata = PredictedFR(fitted.xdata,fitted.kernel,fitted.baseline);
[fitted] = CalculateResidual(fitted);

end