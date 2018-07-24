function [peakiness_mat] = find_pulse_detectors(dff_data_mat_f, stim_mat, stim_mat_simple, sig_cell_mat, frame_time, color_vec, suppress_plots)

n_cells = size(dff_data_mat_f, 2);
odor_list = unique(stim_mat_simple(:, 2));
n_odors = length(odor_list);
odor_dur_list = unique(stim_mat_simple(:, 3));
n_od_durs = length(odor_dur_list);

peakiness_mat = zeros(n_cells, n_odors, n_od_durs) + nan;
for odor_n = 1:n_odors
    odor_ni = odor_list(odor_n);        %actual odor number rather than list count

    for dur_n = 1:(n_od_durs - 1)
        curr_dur = odor_dur_list(dur_n);
        curr_trs = find(stim_mat_simple(:, 2) == odor_ni & stim_mat_simple(:, 3) == curr_dur);
        stim_frs = compute_stim_frs(stim_mat, curr_trs(1), frame_time);

        %calculating auto-correlation to look for pulse responsiveness
        for cell_n = 1:size(dff_data_mat_f, 2)
            n_acq_frs = ceil((stim_mat_simple(curr_trs(1), 7) + (curr_dur.*2.*stim_mat_simple(curr_trs(1), 5)) + stim_mat_simple(curr_trs(1), 10) )./frame_time) - round(5./frame_time);
            curr_ave_trace = mean(dff_data_mat_f(stim_frs(1, 1):(stim_frs(size(stim_frs, 1), size(stim_frs, 2)) + round(5./frame_time)), cell_n, curr_trs), 3, 'omitnan');  %analysing trace only around stimulus period

             %condition to skip cells with insignificant dF/F responses
            if sig_cell_mat(cell_n, odor_ni, dur_n) == 0
                continue
            else
            end

            curr_ave_trace_n = curr_ave_trace - mean(curr_ave_trace);
            curr_dur_fr_num = round(curr_dur./frame_time);
            [autocor, lags] = xcorr(curr_ave_trace_n, curr_dur_fr_num.*10, 'coeff');

            [pks, pksi, w, p] = findpeaks(autocor);

            if suppress_plots == 0
                figure(1)
                plot(lags, autocor)
                hold on
                plot( (pksi  - length(autocor)./2), pks, 'rO')
                hold off
            else
            end

            %getting a measure of peakiness of autocorr at pulse on frequency
            expected_ipi = curr_dur.*2;    %because duration is just the on part of each pulse

            %looping through a range of peak prominence scores to see where pulse timing information disappears from the autocorrelogram
            peakiness_score = nan;
            for p_cutoff = 0:0.1:2            %maximum possible value of prominence is 2, given the maximum range of autocor
                pks(p < p_cutoff) = [];
                pksi(p < p_cutoff) = [];
                w(p < p_cutoff) = [];
                p(p < p_cutoff) = [];
                %skipping if there aren't enough peaks left
                if length(pks) < 3
                    continue
                else
                end                        


                ipi = mean(diff(pksi)).*frame_time;

                %allowing for a 30% error in reporting inter-pulse-interval
                if ipi./expected_ipi > 0.7 && ipi./expected_ipi < 1.3
                    peakiness_score = p_cutoff;
                elseif ipi./(expected_ipi.*2) > 0.7 && ipi./(expected_ipi.*2) < 1.3
                    peakiness_score = p_cutoff;
                else
                end

            end
            peakiness_mat(cell_n, odor_n, dur_n) = peakiness_score;


            if suppress_plots == 0
                figure(2)
                plot(curr_ave_trace)
                stim_frs_i = stim_frs - stim_frs(1, 1);
                add_stim_bar(2, stim_frs_i, color_vec(odor_n, :));
                keyboard
                close figure 2
            else
            end

        end

        if suppress_plots == 0
            
            curr_sig_cells = find(peakiness_mat(:, odor_ni, dur_n) > 0.2);
            curr_insig_cells = 1:1:n_cells;
            curr_insig_cells(curr_sig_cells) = [];
            ave_mat = mean(dff_data_mat_f(:, curr_sig_cells, curr_trs), 3, 'omitnan');

            %sorting by peakiness score
            ave_mat = [curr_sig_cells, ave_mat'];
            ave_mat = sortrows(ave_mat);
            ave_mat = ave_mat(:, 2:end)';

            ave_mat_insig = mean(dff_data_mat_f(:, curr_insig_cells, curr_trs), 3, 'omitnan');

            %plotting
            figure(3)
            imagesc(ave_mat', [0, 1])
            colormap(greymap)
            add_stim_bar(1, stim_frs, color_vec(odor_n, :));
            title('peaky responders')

            figure(4)
            imagesc(ave_mat_insig', [0, 1])
            colormap(greymap)
            add_stim_bar(2, stim_frs, color_vec(odor_n, :));
            title('non-peaky responders')


            close figure 1
            close figure 2
        else
        end


    end
end

