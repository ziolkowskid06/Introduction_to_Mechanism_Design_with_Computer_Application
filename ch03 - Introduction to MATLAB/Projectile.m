%% animate the motion of a projectile being thrown in the air at a 45° angle.

clear variables; close all; clc
 
t = 0:0.1:15;  % [s] vector of times
vx0 = 50;      % [m/s] initial velocity in the x direction
vy0 = 50;      % [m/s] initial velocity in the y direction
g = 9.81;      % [m/s^2] acceleration of gravity

% Open a new plot
figure
% Set a location of a figure:
 % 50 pixels from the bottom of the screen 
 % 50 pixels from the left side of the screen
 % 1200 pixels width and 500 pixels height
set(gcf,'Position',[50 50 1200 500])

% Animation of trajectory
for i = 1:length(t)
  x = vx0 * t(i);                % x position
  y = vy0 * t(i) - 0.5*g*t(i)^2; % y position
 
  plot(x,y,'o','MarkerFaceColor','r')
 
  axis equal
  xlabel('x (m)'); xlim([0 600])
  ylabel('y (m)'); ylim([0 150])
  grid on
  
  % See updates on the screen immediately
  drawnow
end