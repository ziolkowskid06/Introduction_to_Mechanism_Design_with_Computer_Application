%% Makes a plot of a sine wave from zero to 360 degrees.
 
clear variables; close all; clc
 
x = 0:360;
theta = x*pi/180;
 
y = sin(theta);
 
plot(x,y,'LineWidth',2)
xlabel('theta (degrees)'); xticks(0:60:360); xlim([0 360])
ylabel('sin(theta)')
title('The sine function')
grid on

% Save the plot.
saveas(gcf, 'SinePlot.png')