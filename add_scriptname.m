function [] = add_scriptname(fig_n, script_name)
%Syntax: [] = add_scriptname(fig_n, script_name)
%This function inserts a line of text with script_name and the date at the
%bottom of figure fig_n.

%checking for x_label
ax = gca;
xlabel_str = ax.XLabel.String;

figure(fig_n)
inset_vals = get(gca, 'TightInset');

%adding hack for weird inset vals when using imagesc plots
if inset_vals(1, 3) == 0
    inset_vals(1, 3) = 0.03;
    inset_vals(1, 4) = 0.0425;
    imagesc_fig = 1;
else
    imagesc_fig = 0;
end

%figuring out where to put text
ax_vals = axis;
if isempty(xlabel_str) == 1
    x_pos = ax_vals(1, 1);
    
    %accounting for imagesc plot y-axis flip
    if imagesc_fig == 0
        y_pos = ax_vals(1, 3) - 2.*((ax_vals(1, 4) - ax_vals(1, 3))./5);
    elseif imagesc_fig == 1
        y_pos = ax_vals(1, 4) + 1.2.*((ax_vals(1, 4) - ax_vals(1, 3))./5);
    else
    end
    
    %making room for text
    set(gca, 'LooseInset', [inset_vals(1, 1)*1.5, inset_vals(1, 2).*2, inset_vals(1, 3)*1.8, inset_vals(1, 4).*1.5]);

elseif isempty(xlabel_str) == 0
    x_pos = ax_vals(1, 1);
    
    %accounting for imagesc plot y-axis flip
    if imagesc_fig == 0
        y_pos = ax_vals(1, 3) - 2.*((ax_vals(1, 4) - ax_vals(1, 3))./5);
    elseif imagesc_fig == 1
        y_pos = ax_vals(1, 4) + 1.8.*((ax_vals(1, 4) - ax_vals(1, 3))./5);
    else
    end
    
    %making room for text
    set(gca, 'LooseInset', [inset_vals(1, 1).*1.2, inset_vals(1, 2).*1.5, inset_vals(1, 3).*1.5, inset_vals(1, 4).*1.2]);
end

script_name = [script_name, '.m,  ' date];
text(x_pos, y_pos, script_name);