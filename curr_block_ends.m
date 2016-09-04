function [block_ends, block_trs] = curr_block_ends(prot_switch_trials, block_n)
    %This function's output is the first and last trial number of trial
    %block number block_n, as specified by the block ending trials listed
    %in prot_switch_trials.
    %syntax: [block_ends] = curr_block_ends(prot_switch_trials, block_n)
    
    if block_n == 1
        block_ends = [1, prot_switch_trials(block_n)];
        
    elseif block_n > 1
        
        block_ends = [(prot_switch_trials(block_n - 1) + 1), prot_switch_trials(block_n)];
    else
    end
    
    block_trs = block_ends(1):1:block_ends(2);
    
end