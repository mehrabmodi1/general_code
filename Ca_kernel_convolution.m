clear all
close all

%creating a calcium kernel
%CaNew = (1 - fr_time./tau).*Ca + (fr_time./tau)Cbaseline + sp_ht.*n_sp;
%%from Vogelstein, Paninski et al, J neurophys; DOI: 10.1152/jn.01073.2009

tau = .5;           %in s
fr_time = .01;       %in s
Cbaseline = .1;
Ca = Cbaseline;
n_frames = 1000;
sp_ht = 1;
saved_Ca = zeros(1, n_frames);
for t_step = 1:n_frames
    if t_step == 50
        n_sp = 1;
    elseif t_step == 75
        n_sp = .5;
    else
        n_sp = 0;
    end
        
    CaNew = ( (1 - fr_time./tau).*Ca) + ( (fr_time./tau).*Cbaseline) + (sp_ht.*n_sp);
    saved_Ca(1, t_step) = CaNew;
    Ca = CaNew;
    
end

plot(saved_Ca);