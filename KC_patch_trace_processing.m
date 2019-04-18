

TH_path = 'D:\Data\Janelia\Patch\Data_MM\thacq_files\20160802_cell3.mat';

TH_file = load(TH_path);
del = fieldnames(TH_file);      %list of all the fields in the struc TH_file

n_trials = size(del, 1);
clear del
thresh1_val = 3;
thresh2_val = 9.1;
thresh3_val = 0.019;

%Manually identifying thresholds to use for this cell
t_vec = [31:34, (n_trials - 3):n_trials];
for trial_ni = 1:6
    trial_n = t_vec(trial_ni);
    V_trace = eval(['TH_file.Data_' int2str(trial_n) '.data.voltage']).*10;        %in mV because scale of V output from amplifier is 100 mV/mV
    V_trace = double(V_trace);
    T_trace = eval(['TH_file.Data_' int2str(trial_n) '.data.time']);
    %running GUI to manually pick threshold values
    thresh_out = thresh_adjust_gui(V_trace, thresh1_val, thresh2_val, thresh3_val);
    
    thresh1_val = thresh_out(1);
    thresh2_val = thresh_out(2);
    thresh3_val = thresh_out(3);

end


for trial_n = 31:n_trials
    V_trace = eval(['TH_file.Data_' int2str(trial_n) '.data.voltage']).*10;        %in mV because scale of V output from amplifier is 100 mV/mV
    V_trace = double(V_trace);
    T_trace = eval(['TH_file.Data_' int2str(trial_n) '.data.time']);
    
    %running spike detection routine with manually picked threshold values
    spike_s_times = Glenn_DetectSpikes11_m_20160811(V_trace, 100, 1000, thresh1_val, thresh2_val, 1, thresh3_val);
        
    
    
    figure(1)
    plot(V_trace)
%     figure(2)
%     plot(dV2)
%     
    
    
    y_val = max(V_trace).*0.9;
    figure(1)
    hold on
    pks = spike_s_times;
    y_vec = zeros(size(pks, 1), 1) + y_val;
    plot(pks, y_vec, 'r*')
    hold off
    keyboard
      
    
end