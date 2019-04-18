function SpikeSampleTimes = Glenn_DetectSpikes11_m_20160811(V, BPLow, BPHi, NegThresh, PosThresh, PlotFig, diff1_thresh) 
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
time_vec = [1/sf : 1/sf : length(V)/sf] ; % TIME VECTOR FOR PLOTTING

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% FILTER VOLTAGE SIGNAL
% BANDPASS FILTER TO GET RID OF BOTH EXCESSIVELY FAST AND SLOW PARTS OF SIGNAL
filtTrace = bandpassmu(V, sf, BPLow, BPHi) ;
% BOXCAR TO SMOOTH FURTHER
boxc = rectwin(10)/10 ;
filtV = filtfilt(boxc, 1, filtTrace) ;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% TAKE FIRST AND SECOND DERIVATIVES OF FILTERED Vm
diff1 = [0; diff(filtV)] ;
diff2 = [0; diff(diff1)] ;
diff2 = diff2 *-1 ;

% DIPS IN FIRST DERIVATIVE SEEM TO BE MOST DIAGNOSTIC OF SPIKES
% DIPS MUST CROSS CERTAIN # SDs of DERIVATIVE MEASURED IN FIRST 450 msec
% MEASURE VARIABILITY IN FIRST 450 msec OF DERIVATIVE
NegThresh = NegThresh*std(diff1(1:(.45.*sf))) ;
% IDENTIFY CROSSINGS OF THAT SD-BASED THRESHOLD
rr = 1:length(diff1)-1 ;
NegCrossings = find (diff1(rr) > NegThresh & diff1(rr+1) < NegThresh) ;

Candidates = [] ;
tmat = [] ;
PosThresh = PosThresh*std(diff1(1:(.45.*sf))) ;
spike_win = (4./1000).*sf;      % window within which to look for spikes - 4 ms wide
for i = 1:length(NegCrossings)
    if NegCrossings(i) > spike_win + 1
        starts = NegCrossings(i)-spike_win ;
    else starts = 1 ;
    end
        temp = diff1(starts:starts+spike_win) ;
        tmat = [tmat temp] ;
        mtemp = max(temp) ;
        if mtemp > PosThresh
            Candidates = [Candidates NegCrossings(i)];
        end
end

if PlotFig == 1
    figure(3)
    close(gcf)
    hfig = figure(3);
    set(hfig, 'units', 'normalized', 'Position',[0.1 0.1 0.75 0.75]);
    plot(time_vec, diff1, 'b.')
    hold on

    plot((NegCrossings./sf), diff1(NegCrossings), 'rO')
else
end

%%
    if PlotFig == 1
        figure(1); plot(tmat)
    else
    end
    % FIND THE SPIKE PEAKS AS POINTS WHERE 1ST DERIVATIVE SWITCHES FROM POS TO
    % NEG VALUES I.E. PEAKS IN THE INVERTED 2ND DERIVATIVE
    %thresh = 1e-3 ;
    thresh = diff1_thresh;
    SpikeSampleTimes = [] ;
    if isempty(Candidates)
        disp('NoSpikes; stage 1')
    else
        for i = 1:length(Candidates)
            starts = Candidates(i) ;
            %     GRAB THE FIRST DERIVATIVE FOR THE 2msec BEFORE AND AFTER EACH CROSSING
            if starts<2e-3*sf
                %             disp('must be fake spike') % BECAUSE IT COMES BEFORE 200 msec
                %             INTO RECORDING
            else
                temp = diff1(starts-( (2/1000)*sf) : starts + ( (2/1000)*sf) ) ;
                tt = 1:(length(temp) - 1) ;
                %     FIND THE INDEX WHERE diff1 CROSSES ZERO FROM POS TO NEG VALUES
                %     AND MAKE SURE IT RAPIDLY CROSSES ZERO (TO MAKE SURE IT'S REAL)
                try
%                     idx = find (temp(tt) > 0 & temp(tt+1) < 0 ...
%                     & abs (temp(tt)-temp(tt+1)) > thresh ) ;   % [- +]
                    idx = find(abs (temp(tt)-temp(tt+1)) > thresh );

                catch
                    keyboard
                end
                
                    
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
        SpikeTimes = time_vec(1) ;
    else
        SpikeTimes = time_vec(SpikeSampleTimes) ;
    end
    
    %%%%%%%%%%%%%%%
    PlotSoftThresh  = ones(length(time_vec),1) ;
    PlotSoftThresh = PlotSoftThresh * NegThresh ;
    
    PlotHardThresh  = ones(length(time_vec),1) ;
    PlotHardThresh = PlotHardThresh * PosThresh ;
    
    plot_vec_diff_thresh = ones(length(time_vec),1) ;
    plot_vec_diff_thresh = plot_vec_diff_thresh.*diff1_thresh;
    
    % PlotDiff1 = diff1 * 1/range(diff1) ;
    % PlotDiff2 = diff2 * 1/range(diff2) ;
    if PlotFig == 1 ;
        figure(2)
        close(gcf)
        hfig = figure(2);
        set(hfig, 'units', 'normalized', 'Position',[0.1 0.1 0.75 0.75]) ;
        subplot(2,1,1)
        plot(time_vec, V, 'b')
        PlotDiff2 = diff2 * 1/range(diff2) ;
        plot(time_vec, diff1,'b.') ; hold on ;
        plot(time_vec, PlotSoftThresh,'b-') ; hold on ;
        %plot(time_vec, diff2, 'r.') ; hold on ;
        plot(time_vec, PlotHardThresh,'r-') ;
        plot(time_vec, plot_vec_diff_thresh, 'g-');
        
        %         xlim([0.5 .6])
        
        subplot(2,1,2)
        plot(time_vec,filtV-4,'k'); hold on ;
        plot(time_vec,V,'b'); hold on ;
        try
            if size(SpikeSampleTimes, 1) > 0
                plot(SpikeTimes,V(SpikeSampleTimes),'ko', 'markersize', 16) ;
            else
            end
        catch
            keyboard
        end
        %         xlim([0.5 .6])
    end
    
end