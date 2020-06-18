function [train] = rand_train_generator(train_dur, mean_dur, min_pulse_dur, max_pulse_dur)
%This function generates a train of pulses as defined by wait before on and
%wait before off times in s. Each pulse has 
%

% %test vars
% train_dur = 120;
% min_pulse_dur = 0.1;
% max_pulse_dur = 10;
% mean_dur = 2;          %in s, the mean pulse duration

off_dur = 6;        %mean inter-rpulse-interval

%building a reference distribution of pulse durations to draw from 
prob_dist = exprnd(mean_dur, 10000, 1);  %un-sorted, exp distributed random pulse-durs
del = find(prob_dist < min_pulse_dur);
prob_dist(del) = [];                    
del = find(prob_dist > max_pulse_dur);
prob_dist(del) = [];                    


%building a reference distribution of pulse durations to draw from for off
%pulses
prob_dist = exprnd(off_dur, 10000, 1);  %un-sorted, exp distributed random pulse-durs
del = find(prob_dist < min_pulse_dur);
prob_dist(del) = [];                    
del = find(prob_dist > max_pulse_dur);
prob_dist(del) = [];               

%hist(prob_dist)

%building train by drawing pulse-durs from prob_dist

tot_dur = 0;
train = [];
pulse_counter = 1;

while tot_dur < train_dur & abs(train_dur - tot_dur) > 0.001

    if pulse_counter == 1
        %generating the first pulse
        wait_before_on = 0;                 %in s
        rand_n = rand_draw(length(prob_dist));
        wait_before_off = prob_dist(rand_n);
    else
        %generating subsequent pulses
        rand_n = rand_draw(length(prob_dist));
        wait_before_on = prob_dist(rand_n);
        rand_n = rand_draw(length(prob_dist));
        wait_before_off = prob_dist(rand_n);
    end

    train = [train; [wait_before_on, wait_before_off]];
    tot_dur = sum(sum(train));                              %total duration of pulse train so far

    %extending or pruning the last pulse if tot_dur is close to train_dur
    if (train_dur - tot_dur) > 0 && (train_dur - tot_dur) < max_pulse_dur
        t_diff = train_dur - tot_dur;
        train(pulse_counter, 2) = train(pulse_counter, 2) + t_diff;
    elseif (train_dur - tot_dur) < 0
        t_diff = train_dur - tot_dur;
        if train(pulse_counter, 2) > abs(t_diff)
            train(pulse_counter, 2) = train(pulse_counter, 2) + t_diff;
        else
            train(pulse_counter, :) = [];
            pulse_counter = pulse_counter - 1;
        end
    else
    end
    tot_dur = sum(sum(train));

    pulse_counter = pulse_counter + 1;
end

% test plotting    
% 
%t_vec = zeros(train_dur.*1000, 1);
%  t_point = 1;
%  for pulse_n = 1:(pulse_counter - 1)
%      curr_pulse = train(pulse_n, :).*1000;
%      t_vec((t_point + curr_pulse(1)):(t_point + curr_pulse(1) + curr_pulse(2)), 1) = 1;
%      t_point = t_point + curr_pulse(1) + curr_pulse(2);
%      
%  end
% plot(t_vec);
% tot_dur
% figure(2)
% hist(train)

 function [rand_n] = rand_draw(length_vec)
        rand_n = (round(rand(1, 1).*(length_vec - 1)) + 1);
 end
   
end