function [sim_data_mat] = resample_cell_pop(orig_data_mat, n_cells, n_odors, m_sparseness, sd_sparseness, cooperativity)

%BUILDING SIMULATED SIG RESP MAT
sp_vec = normrnd(m_sparseness, sd_sparseness, n_odors, 1);      %vector of randomly generated sparsenesses for each odor with given norm dist props
sp_vec(sp_vec<0) = 0;
n_resp_cells = round(sp_vec.*n_cells);                          %number of sig resp cells for each odor

sig_resp_mat = zeros(n_cells, n_odors);
c_type = sign(cooperativity);
if c_type == 1      %cooperative odor resps
    %perfectly cooperative resps
    for odor_n = 1:n_odors
        n_responders = n_resp_cells(odor_n);
        sig_resp_mat(1:n_responders, odor_n) = 1;
    end
elseif c_type == -1 %anti-cooperative odor resps
    %maximally distributed resps
    for odor_n = 1:n_odors
        n_responders = n_resp_cells(odor_n);
        for responder_n = 1:n_responders
            n_resp_odors_vec = sum(sig_resp_mat, 2);        %vector of the number of odors responded to by each cell
            min_n = min(n_resp_odors_vec);
            min_i = find(n_resp_odors_vec == min_n);        %cell numbers with minimum number of response odors.
            sig_resp_mat(min_i(1), odor_n) = 1;                %picking first minimally responsive cell and making it respond to the current odor.
        end
    end
    
elseif c_type == 0  %randomly distributed odor resps
    for odor_n = 1:n_odors
        n_responders = n_resp_cells(odor_n);
        r_cell_vec = round(rand(n_responders, 1).*(n_cells - 1)) + 1;
        sig_resp_mat(r_cell_vec, odor_n) = 1;
    end
    
else
end

%randomly shuffling some of the initialised odor resps depending on the analog value of cooperativity
if cooperativity ~= 0
    for odor_n = 1:n_odors
        n_responders = n_resp_cells(odor_n);
        n_to_shuff = round(n_responders.*(1 - abs(cooperativity)));            %the number of responder cell numbers picked out of the idealised distributions and randomly re-assigned
        if n_to_shuff <= 0
            continue
        else
        end
        cells_to_shuff = round(rand(n_to_shuff, 1).*(n_responders - 1)) + 1; 
        responder_cells = find(sig_resp_mat(:, odor_n) == 1);
        cells_to_shuff = responder_cells(cells_to_shuff);                   %responder cell numbers picked our randomly whose responses to odor_n will be randomly reassigned
        sig_resp_mat(cells_to_shuff, odor_n) = 0;
        new_rand_resps = round(rand(n_to_shuff, 1).*(n_cells - 1)) + 1;
        sig_resp_mat(new_rand_resps, odor_n) = 1;                               %new, randomly assigned responders to odor_n
    end
else
end
%figures to test initialisation of sig_resp_mat in different cooperativity modes
% figure(1)
% imagesc(sig_resp_mat)
% xlabel('odor number')
% ylabel('cell number')
% 
% figure(2)
% plot(sum(sig_resp_mat, 2))
% xlabel('cell number')
% ylabel('n responded odors')

%%
%Re-sampling real data traces corressponding to sig_resp_mat
n_cells_orig = size(orig_data_mat, 2);
n_frames_orig = size(orig_data_mat(1).traces, 1);
n_reps_orig = size(orig_data_mat(1).traces, 2);
n_odors_orig = size(orig_data_mat(1).traces, 3);
n_durs_orig = size(orig_data_mat(1).traces, 4);

for cell_n = 1:n_cells
    sim_data_mat(cell_n).traces = zeros(n_frames_orig, n_reps_orig, n_odors, n_durs_orig) + nan;
    for odor_n = 1:n_odors
        curr_type = sig_resp_mat(cell_n, odor_n);       %checking if curr cell-odor pair is a responder or not
        
        %randomly picking original data cells until a matching type is
        %found to sample traces from
        sampled = 0;
        while sampled == 0
            curr_cell_orig = round(rand(1, 1).*(n_cells_orig-1)) + 1;
            curr_odor_orig = round(rand(1, 1).*(n_odors_orig-1)) + 1;
            sig_mat_orig = orig_data_mat(curr_cell_orig).info.sig_resps;
            
            curr_type_orig = sum(sig_mat_orig(curr_odor_orig, :));
            if sign(curr_type_orig) == curr_type
                sim_data_mat(cell_n).traces(1:n_frames_orig, 1:n_reps_orig, odor_n, 1:n_durs_orig) = orig_data_mat(curr_cell_orig).traces(:, :, curr_odor_orig, :);
                sampled = 1;
                
            else
            end
            
        end
    end
end
