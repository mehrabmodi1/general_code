function [] = playStack_specint(stack, pause_dur, int_range)
%Syntax: [] = playStack(stack, pause_dur, int_multiplier)
%This function plays back each frame in a 3D stack with a user-defined
%pause (in seconds) between frames and a user-defined colormap scaling factor.

figure(99)
plot_big_fig(99);
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
colormap('gray')
done = 0;
while done == 0
    for frame_n = 1:size(stack, 3)
        frame = squeeze(stack(:, :, frame_n));
        if isempty(int_range) == 0
            curr_int_range = int_range(frame_n, :);
        else
            high_val = 2.*median(reshape(frame, 1, []));
            curr_int_range = [0, high_val];
        end
        
        high_int = curr_int_range(1, 2);
        if high_int == 0 || isnan(high_int) == 1
            high_int = 1;
        else
        end
        try
            imagesc(frame, curr_int_range);
        catch
            keyboard
        end
        title(['Frame ' int2str(frame_n), ' of  ', int2str(size(stack, 3)), '.'])
        drawnow
        pause(pause_dur.*0.001)
    end
    choice = questdlg('Play again?', 'Playing stack', 'Re-play', 'Stop', 'Stop');
    if strcmp(choice, 'Stop') == 1
        done = 1;
    else
    end
    
end
close figure 99