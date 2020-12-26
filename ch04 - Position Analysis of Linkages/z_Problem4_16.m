% Problem 4.16 (inverted slider-crank linkage)
% by Damian Ziolkowski, December 22, 2020
 
% Prepare Workspace
clear variables; close all; clc;
 
% Linkage dimensions
a = 0.030;           % crank length (m)
c = 0.080;           % rocker length (m)
d = 0.090;           % length between ground pins (m)
p = 0.200;           % slider length (m)
e = 0.050;           % PQ length (m)
delta = 45*pi/180;   % angle between slider and rocker (converted to rad)
fix = -90*pi/180;    % fixed CPQ angle
eps = 0.000001;      % tiny number added to theta2 to keep it in bounds
h = c*sin(delta);    % h is a constant, only calculate it once
 
% Ground pins
x0 = [0;0];   % ground pin at A (origin)
xD = [d;0];   % ground pin at D

% Grashof check
if ((a + d) < c*sin(delta))
  disp('Linkage cannot be assembled')
  return
else
  if (abs(a - d) >= c*sin(delta))
    disp('Linkage is Grashof')
    theta2min = 0;
    theta2max = 2*pi;
  else
    disp('Linkage is not Grashof')
    theta2min = acos((a^2 + d^2 - (c*sin(delta))^2)/(2*a*d)) + eps;
    theta2max = 2*pi - theta2min;
  end
end

N = 361;   % number of times to perform position calculations
[xB, xC, xP, xQ] = deal(zeros(2, N));          % allocate space for positions
[theta2, theta3, theta4] = deal(zeros(1, N)); % allocate space for angles
                        % "deal" distributes inputs to outputs

% Perform calculations for each angle
for i = 1:N
  theta2(i) = (i-1)*(theta2max - theta2min)/(N-1) + theta2min;
  r = d - a*cos(theta2(i));
  s = a*sin(theta2(i));
  f2 = r^2 + s^2;  % f squared
 
  b(i) = c * cos(delta) + sqrt(f2 - h^2);
  g = b(i) - c*cos(delta);
 
  theta3(i) = atan2((h*r - g*s),(g*r + h*s));
  theta4(i) = theta3(i) + delta;
 
 % Calculate unit vectors
  [e2,n2]   = UnitVector(theta2(i));
  [e3,n3]   = UnitVector(theta3(i));
  [e4,n4]   = UnitVector(theta4(i));
  [ePQ,nPQ] = UnitVector(theta3(i) + fix);
 
 % Solve for positions of points B, C and P on the linkage
  xB(:,i) = FindPos(     x0, a, e2);
  xC(:,i) = FindPos(     xD, c, e4);
  xP(:,i) = FindPos(xB(:,i), p, e3);
  xQ(:,i) = FindPos(xP(:,i), e, ePQ);
end

% Plot paths for points P, and Q
plot(xP(1,:),xP(2,:), 'Color', [204/255 102/255 0/255])
hold on
plot(xQ(1,:),xQ(2,:), 'Color', [0/255 153/255 76/255])
 
% Specify angle at which to plot linkage
iTheta = 120;
 
% Plot linkage
plot([       x0(1) xB(1,iTheta)],[       x0(2) xB(2,iTheta)],...
      'Linewidth',2,'Color','k');
plot([xB(1,iTheta) xP(1,iTheta)],[xB(2,iTheta) xP(2,iTheta)],...
      'Linewidth',2,'Color','k');
plot([xD(1)        xC(1,iTheta)],[xD(2)        xC(2,iTheta)],...
      'Linewidth',2,'Color','k');
plot([xP(1,iTheta) xQ(1,iTheta)],[xP(2,iTheta) xQ(2,iTheta)],...
      'Linewidth',2,'Color','k');
 
% Plot joints on linkage
 plot([x0(1) xD(1) xB(1,iTheta) xC(1,iTheta) xP(1,iTheta) xQ(1,iTheta)],...
      [x0(2) xD(2) xB(2,iTheta) xC(2,iTheta) xP(2,iTheta) xQ(2,iTheta)],...
      'o','MarkerSize',5,'MarkerFaceColor','k','Color','k');
 
% Plot the labels of each joint
text(xB(1,iTheta),xB(2,iTheta)+.010,'B','HorizontalAlignment','center');
text(xC(1,iTheta),xC(2,iTheta)+.010,'C','HorizontalAlignment','center');
text(xP(1,iTheta),xP(2,iTheta)+.010,'P','HorizontalAlignment','center');
text(xQ(1,iTheta),xQ(2,iTheta)-.010,'Q','HorizontalAlignment','center');
 
axis equal
grid on
 
title('Paths of points P and Q on Inverted Slider-Crank')
xlabel('x-position [m]')
ylabel('y-position [m]')
legend('Point P', 'Point Q','Location','SouthEast') 

% Save the plot
saveas(gcf, 'z_Problem4_16 - plot.png')