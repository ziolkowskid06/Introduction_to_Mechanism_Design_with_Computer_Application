% Problem 4.11 (fourbar non-grashof linkage)
% by Damian Ziolkowski, December 19, 2020

% Prepare Workspace
clear variables; close all; clc;

% Linkage dimensions
a = 0.050;         % crank length (m)
b = 0.075;         % coupler length (m)
c = 0.085;         % rocker length (m)
d = 0.115;         % length between ground pins (m)
 
% Ground pins
x0 = [0; 0];   % ground pin at A (origin)
xD = [d; 0];   % ground pin at D

% Grashof Check
S = min([a b c d]);  % length of shortest link
L = max([a b c d]);  % length of longest link
T = sum([a b c d]);  % total of all link lengths
PQ = T - S - L;  % length of P plus length of Q
if (S+L < PQ)  % Grashof condition
  disp('Linkage is Grashof.')
  theta2min = 0;
  theta2max = 2*pi;
else  % if not Grashof, terminate program
  disp('Linkage is not Grashof')
  theta2max = acos((a^2 + d^2 - (b + c)^2)/(2*a*d));
  theta2min = -theta2max;
end

N = 361;   % number of times to perform position calculations
[xB, xC] = deal(zeros(2, N));          % allocate space for points
[theta2, theta3, theta4] = deal(zeros(1, N)); % allocate space for angles
                        % "deal" distributes inputs to outputs

% Perform calculations for each angle
for i = 1:N
  theta2(i) = (i-1)*(theta2max - theta2min)/(N-1) + theta2min;  % crank angle
 
 % Conduct position analysis to solve for theta3 and theta4
  r = d - a*cos(theta2(i));
  s = a*sin(theta2(i));
  f2 = r^2 + s^2;                     % f squared
  delta = acos((b^2+c^2-f2)/(2*b*c)); % angle between coupler and rocker
 
  g = b - c*cos(delta);
  h = c*sin(delta);
 
  theta3(i) = atan2((h*r - g*s),(g*r + h*s));  % coupler angle
  theta4(i) = theta3(i) + delta;               % rocker angle
 
 % Calculate unit vectors
  [e1,n1] = UnitVector(0);
  [e2,n2] = UnitVector(theta2(i));
  [e3,n3] = UnitVector(theta3(i));
  [e4,n4] = UnitVector(theta4(i));
 
 % Solve for positions of points B, C and P on the linkage
  xB(:,i) = FindPos(     x0,   a,   e2);
  xC(:,i) = FindPos(     xD,   c,   e4);
end

% Plot paths for points B, and C
plot(xB(1,:),xB(2,:), 'Color', [0/255 128/255 255/255])
hold on
plot(xC(1,:),xC(2,:), 'Color', [204/255 102/255 0/255])
axis equal
grid on

% Specify angle at which to plot linkage
iTheta = 320;
  
% Plot linkage
plot([x0(1) xB(1,iTheta)],[x0(2) xB(2,iTheta)],'Linewidth',2,'Color','k');
plot([xB(1,iTheta) xC(1,iTheta)],[xB(2,iTheta) xC(2,iTheta)],'Linewidth',2,'Color','k');
plot([xD(1) xC(1,iTheta)],[xD(2) xC(2,iTheta)],'Linewidth',2,'Color','k');
 
% Plot joints on linkage
plot([x0(1) xD(1) xB(1,iTheta) xC(1,iTheta) ],...
     [x0(1) xD(2) xB(2,iTheta) xC(2,iTheta) ],...
     'o','MarkerSize',5,'MarkerFaceColor','k','Color','k');
 
% Plot the labels of each joint
text(       x0(1),       x0(2)-.005,'A','HorizontalAlignment','center');
text(       xD(1),       xD(2)-.005,'D','HorizontalAlignment','center');
text(xB(1,iTheta),xB(2,iTheta)+.005,'B','HorizontalAlignment','center');
text(xB(1,iTheta),xB(2,iTheta)+.005,'B','HorizontalAlignment','center');
text(xC(1,iTheta),xC(2,iTheta)+.005,'C','HorizontalAlignment','center');
axis equal
grid on
title('Paths of points B, and C on Fourbar Non-Grashof Linkage')
xlabel('x-position [m]')
ylabel('y-position [m]')
legend('Point B', 'Point C','Location','SouthEast')

% Save the plot
saveas(gcf, 'z_Problem4_11 - plot.png')