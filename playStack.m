function [] = playStack(stack, pause_dur, int_multiplier)
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
        high_int = double(max(max(frame))).*0.6.*int_multiplier;
        if high_int == 0 || isnan(high_int) == 1
            high_int = 1;
        else
        end
        try
            imagesc(frame, [0, high_int]);
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