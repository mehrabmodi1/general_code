function [is_conn] = check_ser_connected(portname)

p_list = instrfind;

n_devs = size(p_list, 2);
is_conn = 0;
for dev_n = 1:n_devs
    curr_name = p_list(dev_n).Name(8:end);
    
    if strcmp(curr_name, portname) == 1
        is_conn = 1;
        break
    else
    end
    
end
