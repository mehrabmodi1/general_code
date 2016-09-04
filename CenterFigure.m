% Function to center a figure on the screen (if it's not maximized).
function CenterFigure(handleToFigure)
% The figure Position property
% does not include the window borders, so this example uses a width of 5 pixels
% on the sides and bottom and 30 pixels on the top.
borderWidth = 5;
titleBarWidth = 30;
% Ensure root units are pixels and get the size of the screen:
set(0, 'Units', 'pixels');
set(handleToFigure, 'Units', 'pixels');
% Get the screen size in pixels.
screenSize = get(0,'ScreenSize');
% Get the size of the window.
initialFigurePosition = get(handleToFigure, 'Position');
% Create an array that will center it.
centeredX = (screenSize(3) - initialFigurePosition(3)) / 2;
centeredY = (screenSize(4) - initialFigurePosition(4)) / 2;
centeredPosition = [centeredX,...
centeredY,...
initialFigurePosition(3),...
initialFigurePosition(4)];
set(handleToFigure, 'Position', centeredPosition);
return; % from CenterFigure()