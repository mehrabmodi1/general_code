function [] = KCpop_odor_resp_plotter(dff_data_mat_f, stim_mat_simple, stim_mat, sig_cell_mat, frame_time, start_fig_n, script_name, od_col_ns, dur_col_ns)

color_vec = load('C:\Users\Mehrab\Google Drive\Backup\Stuff\CSHL\Glenn lab\Code\std_color_vec.txt');
a = colormap('bone');
greymap = flipud(a);
[del, odor_names] = xlsread('C:\Data\Code\general_code_old\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList.xls', 1);


odor_list = unique(stim_mat_simple(:, od_col_ns(1)));
odor_dur_list = unique(stim_mat_simple(: ,dur_col_ns(1)));
if isempty(sig_cell_mat) == 0
    %identifying all cells that are sig responders to at least one stimulus
    sig_cell_mat = sum(sig_cell_mat, 3);
    sig_cell_mat = sum(sig_cell_mat, 2);
    sig_cells = find(sig_cell_mat > 0);
elseif isempty(sig_cell_mat) == 1
    sig_cells = 1:1:size(dff_data_mat_f, 2);
else
end

for odor_n = 1:length(odor_list)
    odor_ni = odor_list(odor_n);
    for dur_n = 1:length(odor_dur_list)
        od_duri = odor_dur_list(dur_n);
        
        if od_duri == 0.1
            continue
        else
        end
        
        curr_trs = find(stim_mat_simple(:, od_col_ns(1)) == odor_ni & stim_mat_simple(:, dur_col_ns(1)) == od_duri);
        
        %stim_frs = compute_stim_frs(stim_mat, curr_trs(1), frame_time);
        stim_frs = compute_stim_frs_modular(stim_mat, curr_trs(1), frame_time);
        stim_frs = stim_frs{1};    
        curr_traces = mean(dff_data_mat_f(:, sig_cells, curr_trs), 3, 'omitnan');      %ave KC pop resp traces to curr od-dur pair
        curr_od_name = odor_names{odor_ni};
        if exist('max_vals') == 0
            max_vals = max(curr_traces, [], 1);
        else
        end
        
        %sorting traces by pk responses in od1-dur1.
        curr_traces = ([max_vals; curr_traces])';
        curr_traces = sortrows(curr_traces);
        curr_traces(:, 1) = [];
        %curr_traces(31, :) = [];
        figure(start_fig_n)
        imagesc(curr_traces, [0, 1.5])
        colormap(greymap)
        ylabel('cell number')
        set_xlabels_time(start_fig_n, frame_time, 10)
        fig_wrapup_mod(start_fig_n, 'tall', script_name)   
        add_stim_bar(start_fig_n, stim_frs, [0.5, 0.5, 0.5])
                     
        disp([num2str(start_fig_n), ' ', curr_od_name, ', ' num2str(od_duri), 's'])
        start_fig_n = start_fig_n + 1;
    end
end
