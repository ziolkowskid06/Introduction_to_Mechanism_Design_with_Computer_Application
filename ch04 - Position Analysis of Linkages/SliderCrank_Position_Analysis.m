% Performs a position analysis on the slider-crank linkage and
% plots the piston position as a function of crank angle
% by Damian Ziolkowski, November 26, 2020

% Prepare Workspace
clear variables; close all; clc;

% Linkage dimensions
a = 0.040;         % crank length (m)
b = 0.120;         % connecting rod length (m)
c = 0.0;           % vertical slider offset (m)
 
% Ground pins
x0 = [0;0];        % ground pin at A (origin)

N = 361;   % number of times to perform position calculations
[xB,xC] = deal(zeros(2,N));          % allocate space for pins B and C
[theta2,theta3,d] = deal(zeros(1,N));% allocate space for link angles
                 % "deal" distributes inputs to outputs

% Perform calculations for each angle
for i = 1:N
  theta2(i) = (i-1)*(2*pi)/(N-1);
  theta3(i) = asin((c - a*sin(theta2(i)))/b);
  d(i) = a*cos(theta2(i)) + b*cos(theta3(i));
 
% Calculate unit vectors
  [e1,n1] = UnitVector(0);
  [e2,n2] = UnitVector(theta2(i));
  [e3,n3] = UnitVector(theta3(i));
  
% Solve for position of point B on the linkage
  xB(:,i) = FindPos(     x0,   a,   e2);
% Not needed but validate the position of point C
  xC(:,i) = FindPos(xB(:,i),   b,   e3);  
end

% Plot the piston position
plot(theta2*180/pi, d*1000, 'o', 'Color', [255/255 155/255 155/255])
hold on
plot(theta2*180/pi, xC(1,:)*1000, 'Color', [0/255 153/255 76/255])
title('Piston Position versus Crank Angle for Slider-Crank')
xlabel('Crank angle (degrees)')
ylabel('Position (mm)')
legend('d','xC')
grid on
set(gca,'xtick',0:60:360)
xlim([0 360])

% Save the plot
saveas(gcf, 'SliderCrank_Position_Analysis - plot.png')