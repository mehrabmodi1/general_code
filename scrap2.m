clear all
close all
figure(1)

plot(rand(10, 10).*100)
script_name = mfilename;
xlabel('xlabel')
fig_wrapup(1, script_name);