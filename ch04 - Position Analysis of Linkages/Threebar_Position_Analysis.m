% Conducts a position analysis on the threebar crank-slider linkage
% by Damian Ziolkowski, November 25, 2020

% Prepare Workspace
clear variables; close all; clc;
 
% Linkage dimensions
a = 0.100;         % crank length (m)
d = 0.150;         % length between ground pins (m)
p = 0.300;         % slider length (m)
 
% Ground pins
xA = [0;0];        % point A (the origin)
xD = [d;0];        % point D
 
N = 361;   % number of times to perform position calculations
[xB,xP]           = deal(zeros(2, N)); % allocate space for position of B,P
[theta2,theta3,b] = deal(zeros(1, N)); % allocate space for link angles
                 % "deal" distributes inputs to outputs
 
% Perform calculation for each angle
for i = 1:N
  % Map iterations to angles to make sure crank makes complete rotation  
  theta2(i) = (i-1)*(2*pi)/(N-1);
  % Use "atan2" to calculate angle for all four quadrants
  theta3(i) = atan2(-a*sin(theta2(i)), d - a*cos(theta2(i)));
  b(i) = (d - a*cos(theta2(i)))/cos(theta3(i));
 
% Calculate unit vectors
  [e2,n2]   = UnitVector(theta2(i));
  [e3,n3]   = UnitVector(theta3(i));
 
% Solve for positions of points B and P on the linkage
  xB(:,i) = FindPos(     xA, a, e2);
  xP(:,i) = FindPos(xB(:,i), p, e3);
end

% Path of point B
plot(xB(1,:), xB(2,:), 'Color', [153/255 153/255 153/255])
hold on
% Path of point P
plot(xP(1,:), xP(2,:), 'Color', [0/255 153/255 76/255])
 
% Specify angle at which to plot a linkage
iTheta = 80;
 
% Plot crank and slider
plot([xA(1)        xB(1,iTheta)],...
     [xA(2)        xB(2,iTheta)], 'Linewidth', 2, 'Color', 'k');
plot([xB(1,iTheta) xP(1,iTheta)],...
     [xB(2,iTheta) xP(2,iTheta)], 'Linewidth', 2, 'Color', 'k');
 
% Plot joints on linkage
plot([xA(1) xD(1) xB(1,iTheta) xP(1,iTheta)],...
     [xA(2) xD(2) xB(2,iTheta) xP(2,iTheta)],...
     'o', 'MarkerSize', 5, 'MarkerFaceColor', 'k', 'Color', 'k');
 
% Plot the labels of each joint
text( xA(1)-0.015,              xA(2), 'A', 'HorizontalAlignment', 'center');
text(xB(1,iTheta), xB(2,iTheta)+0.015, 'B', 'HorizontalAlignment', 'center');
text(       xD(1),        xD(2)+0.015, 'D', 'HorizontalAlignment', 'center');
text(xP(1,iTheta), xP(2,iTheta)+0.015, 'P', 'HorizontalAlignment', 'center');
title('Paths of points B and P on the Threebar Linkage')
xlabel('x-position [m]')
ylabel('y-position [m]')
legend('Point B', 'Point P', 'Location', 'SouthEast')
axis equal
grid on 

% Save the plot
saveas(gcf, 'Threebar_Position_Analysis - plot.png')