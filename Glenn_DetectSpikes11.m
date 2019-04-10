function SpikeSampleTimes = DetectSpikes10(V, BPLow, BPHi, NegThresh, PosThresh, PlotFig) ;
% FILTERS VOLTAGE TRACE AND USES NEGATIVE PEAKS IN FIRST DERIVATIVE TO FIND
% SPIKES SINCE THE REPOLARIZATION OF THE SPIKE IS WHAT IS MOST DIAGNOSTIC

% TO IMPROVE: NORMALIZE THE VM AND THE DERIVATIVES AND PLOT ON TOP OF ONE ANOTHER
% NEED TO USE SECOND DERIVATIVE TO FIND SPIKES IN VOLTAGE TO ACCURATELY GET SPIKES
% IMPLEMENT A BIG FIRST DERIVATIVE PEAK THAT PRECEEDS THE FIRST DERIVATIVE
% DIP - EITHER THAT OR INCREASE THE THRESHOLD FOR THE SPEED OF THE ZERO CROSSING
% WOULD BE HELPFUL TO PLOT THE MAGNITUDES OF FIRST DERIV PEAKS AND DIPS TO SEE WHERE TO DRAW THE THRESHOLD
% MIGHT ALSO BE USEFUL TO SCAN THE DERIVATIVES AGAINST A TEMPATE.

% INITIALIZE VARIABLES
sf = 10000 ;  % SAMPLING FREQUENCY
time = [1/sf : 1/sf : length(V)/sf] ; % TIME VECTOR FOR PLOTTING

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% FILTER VOLTAGE SIGNAL
% BANDPASS FILTER TO GET RID OF BOTH EXCESSIVELY FAST AND SLOW PARTS OF SIGNAL
filtTrace = bandpassmu(V, sf, BPLow, BPHi) ;
% BOXCAR TO SMOOTH FURTHER
boxc = boxcar(10)/10 ;
filtV = filtfilt(boxc, 1, filtTrace) ;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% TAKE FIRST AND SECOND DERIVATIVES OF FILTERED Vm
diff1 = [0; diff(filtV)] ;
diff2 = [0; diff(diff1)] ;
diff2 = diff2 *-1 ;

% DIPS IN FIRST DERIVATIVE SEEM TO BE MOST DIAGNOSTIC OF SPIKES
% DIPS MUST CROSS CERTAIN # SDs of DERIVATIVE MEASURED IN FIRST 450 msec
% MEASURE VARIABILITY IN FIRST 450 msec OF DERIVATIVE
NegThresh = NegThresh*std(diff1(1:4500)) ;
% IDENTIFY CROSSINGS OF THAT SD-BASED THRESHOLD
rr = 1:length(diff1)-1 ;
NegCrossings = find (diff1(rr) > NegThresh & diff1(rr+1) < NegThresh) ;

Candidates = [] ;
tmat = [] ;
PosThresh = PosThresh*std(diff1(1:4500)) ;
for i = 1:length(NegCrossings)
    if NegCrossings(i) >41
        starts = NegCrossings(i)-40 ;
    else starts = 1 ;
    end
        temp = diff1(starts:starts+40) ;
        tmat = [tmat temp] ;
        mtemp = max(temp) ;
        if mtemp > PosThresh
            Candidates = [Candidates NegCrossings(i)]
        end
    end
  figure;plot(tmat)  
    % FIND THE SPIKE PEAKS AS POINTS WHERE 1ST DERIVATIVE SWITCHES FROM POS TO
    % NEG VALUES I.E. PEAKS IN THE INVERTED 2ND DERIVATIVE
    thresh = 1e-3 ;
    SpikeSampleTimes = [] ;
    if isempty(Candidates)
        disp('NoSp')
    else
        for i = 1:length(Candidates)
            starts = Candidates(i) ;
            %     GRAB THE FIRST DERIVATIVE FOR THE 2msec PRIOR TO EACH CROSSING
            if starts<2000
                %             disp('must be fake spike') % BECAUSE IT COMES BEFORE 200 msec
                %             INTO RECORDING
            else
                temp = diff1(starts-(2/1000*sf) : starts) ;
                tt = 1:20 ;
                %     FIND THE INDEX WHERE diff1 CROSSES ZERO FROM POS TO NEG VALUES
                %     AND MAKE SURE IT RAPIDLY CROSSES ZERO (TO MAKE SURE IT'S REAL)
                idx = find (temp(tt) > 0 & temp(tt+1) < 0 ...
                    & abs (temp(tt)-temp(tt+1)) > thresh ) ;   % [- +]
                if isempty(idx)
                    %         disp('No Vpeak detected')
                else
                    % TAKE THE INDEX CLOSEST TO THE THRESHOLD CROSSING.  THIS SHOULD BE ADJUSTED SINCE I WANT TO USE SECOND DERIVATIVE TO FIND VPEAK OF SPIKE
                    idx = idx(end) ;
                    SpikeSampleTimes = [SpikeSampleTimes Candidates(i)+idx] ;
                end
            end
        end
    end
    if isempty(SpikeSampleTimes)
        disp('NoSpikes')
        SpikeTimes = time(1) ;
    else
        SpikeTimes = time(SpikeSampleTimes) ;
    end
    
    %%%%%%%%%%%%%%%
    PlotSoftThresh  = ones(length(time),1) ;
    PlotSoftThresh = PlotSoftThresh * NegThresh ;
    
    PlotHardThresh  = ones(length(time),1) ;
    PlotHardThresh = PlotHardThresh * PosThresh ;
    
    % PlotDiff1 = diff1 * 1/range(diff1) ;
    % PlotDiff2 = diff2 * 1/range(diff2) ;
    if PlotFig == 1 ;
        figure('Position',[0 800 1200 800]) ;
        subplot(2,1,1)
        plot(time, V, 'b')
        PlotDiff2 = diff2 * 1/range(diff2) ;
        plot(time, diff1,'b.') ; hold on ;
        plot(time, PlotSoftThresh,'b-') ; hold on ;
        plot(time, diff2, 'r.') ; hold on ;
        plot(time, PlotHardThresh,'r-') ;
        %         xlim([0.5 .6])
        
        subplot(2,1,2)
        plot(time,filtV-4,'k'); hold on ;
        plot(time,V,'b'); hold on ;
        plot(SpikeTimes,V(SpikeSampleTimes),'ko', 'markersize', 16) ;
        %         xlim([0.5 .6])
    end
end