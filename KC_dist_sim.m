function [] = KC_dist_sim()

%specifying basic params

%obtained odor response fractions from Rob's paper
frac_list = [0.01, 0.01, 0.01, 0.025, 0.01, 0.01, 0.04, 0.04, 0.02, 0.025, 0.12, 0.01, 0.04, 0.05, 0.09, 0.02, 0.04, 0.05, 0.08, 0.03, 0.04, 0.04, 0.11, 0.05, 0.055, 0.07, 0.09, 0.10, 0.11, 0.02, 0.04, 0.05, 0.08, 0.1, 0.1, 0.11, 0.03, 0.06, 0.065, 0.11, 0.12, 0.17];
n_fracs = size(frac_list, 2);
%histogram(nKC_dist)
n_runs = 100;
n_odors = 500;
nKCs = 1000;
mean_wt_high = 0.9;
mean_wt_low = 0.001;
sd_wt_high = 0;       %multiplier of mean_wt
sd_wt_low = 0;        %multiplier of mean_wt
plast_step = 0.85;
sd_wt_plast = 0;      %multiplier of mean_wt
sparseness_multiplier = 10;
make_plots = 1;         %set to 0 to prevent plotting

%running the simulation
%loop to do many repeated draws of odor response distributions
for iteration_n = 1:n_runs
    all_responders = [];
    %Setting up baseline state
    n_responders_list = [];
    for odor_n = 1:n_odors
        %setting up responder KC vector for current odor
        fraci = round(rand(1, 1).*(n_fracs - 1)) + 1;
        n_responders = min([(frac_list(fraci).*nKCs).*sparseness_multiplier], nKCs);
        n_responders_list = [n_responders_list; n_responders];
        
        KC_list = 1:1:nKCs;
        r_list = randperm(nKCs);
        KC_list = [r_list', KC_list'];
        KC_list = sortrows(KC_list);                %this is a randomly re-ordered list of cell indices
        responder_list = KC_list(1:n_responders, 2);
        odor_data(odor_n).KC_list = responder_list;
        odor_data(odor_n).learned = 0;
        
        all_responders = [all_responders; responder_list];
        
        %assigning weights to pot_MBON and dep_MBON
        wt_matrix(1:nKCs, 1) = (sd_wt_high.*mean_wt_high).*randn(nKCs, 1) + mean_wt_high;      %baseline wts for the dep MBON
        wt_matrix(1:nKCs, 2) = (sd_wt_low.*mean_wt_low).*randn(nKCs, 1) + mean_wt_low;        %baseline wts for the pot MBON
        
        odor_data_base = odor_data;
        
        %computing baseline MBON activation and saving to build baseline distributions
        MBON_pot_resps(odor_n, 1) = sum(wt_matrix(odor_data(odor_n).KC_list, 2).*1);
        MBON_dep_resps(odor_n, 1) = sum(wt_matrix(odor_data(odor_n).KC_list, 1).*1);

    end
    
   
    [bin_vec, del, del2] = plot_resp_dists(MBON_pot_resps, MBON_dep_resps, [], make_plots, n_odors);
    
    %loop to simulate learning multiple odors
    learning_order = randperm(n_odors);    
    for l_odor_n = 1:n_odors
        %disp(['learned odor n ' int2str(l_odor_n)])
        odor_data(l_odor_n).learned = 1;
        
        curr_odor = learning_order(l_odor_n);
        curr_KCs = odor_data(l_odor_n).KC_list;
                
        %generating rand wt steps to add to the KCs that represent curr_odor
        high_wts = (sd_wt_plast.*plast_step).*randn(length(curr_KCs), 1) + plast_step;
        low_wts = (sd_wt_plast.*plast_step).*randn(length(curr_KCs), 1) + plast_step;
        
        %assigning new wts
        wt_matrix(curr_KCs, 1) = max([(wt_matrix(curr_KCs, 1) - low_wts), (sd_wt_low.*randn(length(curr_KCs), 1) + mean_wt_low)], [], 2);      %updating wts for the dep MBON
        wt_matrix(curr_KCs, 2) = min([(wt_matrix(curr_KCs, 2) + high_wts), (sd_wt_high.*randn(length(curr_KCs), 1) + mean_wt_high)], [], 2);     %updating wts for the pot MBON
        %keyboard
        %recomputing MBON responses
        MBON_pot_resps(l_odor_n, 1) = sum(wt_matrix(odor_data(l_odor_n).KC_list, 2).*1);
        MBON_dep_resps(l_odor_n, 1) = sum(wt_matrix(odor_data(l_odor_n).KC_list, 1).*1);
    
        plot_resp_dists(MBON_pot_resps, MBON_dep_resps, bin_vec, make_plots, n_odors)
    
    end
    
    
    
    
  keyboard  
end

end

function [bin_vec, dep_his, pot_his] = plot_resp_dists(MBON_pot_resps, MBON_dep_resps, bin_vec, make_plots, n_odors)
    if isempty(bin_vec) == 1
        max_val = max([MBON_pot_resps; MBON_dep_resps]);
        min_val = min([MBON_pot_resps; MBON_dep_resps]);
        bin_vec = linspace(min_val, max_val, 100);
    else
    end
    
    dep_his = hist(MBON_dep_resps, bin_vec);
    pot_his = hist(MBON_pot_resps, bin_vec);
    
    if make_plots == 1
        figure(1)
        plot(bin_vec, dep_his, 'b', 'LineWidth', 2)
        hold on
        plot(bin_vec, pot_his, 'r', 'LineWidth', 2)
        hold off
        xlabel('MBON excitation')
        ylabel('counts')
        axis([0, max(bin_vec), 0, n_odors])
        drawnow
        
        %pause(.1)
    else
    end
end

