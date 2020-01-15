function pairing_tr_n_vec = find_pairing_tr(params)
pairing_tr_n_vec = [];
for tr_n = 1:size(params, 2)
    if params(tr_n).led_on == 1
        pairing_tr_n_vec = [pairing_tr_n_vec; tr_n];
    else
    end
end


