% Problem 4.5 (ant trajectory)
% by Damian Ziolkowski, December 19, 2020

% Prepare Workspace
clear variables; close all; clc;

% Data
t = 10;                 % time (s)
N = 0:0.01:t;           % number of iteration  

% Ground pin
x0 = [0;0];        % ground pin at A (origin)

% Allocate spaces for variables
     [xANT] = deal(zeros(2,length(N)));        % allocate space for point B and Ant coordinates
[theta2, s] = deal(zeros(1,length(N)));        % allocate space for angle and Ant distance

% Main loop
for i = N
   ind = find(N==i); 
   theta2(ind) = (4*i);                  % angle (radians)
   s(ind) = 10*i;                        % ant distance (mm)
   
   % Calculate unit vector
   [eANT, nANT] = UnitVector(theta2(ind));
   
   % Ant coordinates
   xANT(:,ind) = FindPos(x0, s(ind), eANT);  
end

% Plot trajectory for the ant
plot(xANT(1,:),xANT(2,:),'k')
title('Ant trajectory for 10 seconds')
xlabel('x-position [mm]')
ylabel('y-position [mm]')
axis equal
grid on

% Save the plot
saveas(gcf, 'z_Problem4_05 - plot.png')