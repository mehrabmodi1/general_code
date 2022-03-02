curr_dir = 'C:\Data\Data\Analysed_data\Suite2P_results\20200622\fly2_c305a_KC_handover_train\1\';   %A'/B'
curr_dir = 'C:\Data\Data\Analysed_data\Suite2P_results\20200804\fly1_c739XopGC6f\1\';   %A/B

[stim_mat, stim_mat_simple, column_heads, color_vec, good_tr_list, params_orig, PID_traces] = load_params_trains_modular(curr_dir, [], 0.091);    %reading in trial stimulus parameters after matching time stamps to F traces
        