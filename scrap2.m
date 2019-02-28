clear all
close all
figure(1)

imagesc(rand(10, 10).*100)
script_name = mfilename;
xlabel('xlabel')
ylabel('ylabel')
fig_wrapup(1, script_name)