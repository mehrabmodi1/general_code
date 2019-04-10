function data = aq_data_bk(src, event)
%This function is called by the listener in the background acquisition to
%acquire the PID signal and write it to file.

data = [event.Data, event.TimeStamps];
curr_dir = curr_aq_direc;
save([curr_dir '\last_PID_trace.mat'], 'data');
