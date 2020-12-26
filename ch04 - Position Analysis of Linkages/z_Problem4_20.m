% Problem 4.20 (claw mechanism)
% by Damian Ziolkowski, December 25, 2020

% Prepare Workspace
clear variables; close all; clc;

% Linkage dimensions
a = 40;                         % AB crank length (mm)
b = 60;                         % BC coupler length (mm)
c = 60;                         % DC rocker length (mm)
d = 60;                         % AD length between ground pins (mm)
u = 60;                         % CP length (mm)
v = 60;                         % CR length (mm)
y = 20;                         % PQ length (mm)
z = 20;                         % RS length (mm)
gamma2 = -(180-70)*pi/180;      % CPQ angle (rad) (CW rotation)
gamma4 =  (180-70)*pi/180;      % CRS angle (rad) (CCW rotation)

 
% Ground pins
x0 = [0; 0];   % ground pin at A (origin)
xD = [d; 0];   % ground pin at D

% Grashof Check
S = min([a b c d]);  % length of shortest link
L = max([a b c d]);  % length of longest link
T = sum([a b c d]);  % total of all link lengths
PQ = T - S - L;      % length of P plus length of Q
if (S+L < PQ)        % Grashof condition
  disp('Linkage is Grashof.')
  theta2min = 0;
  theta2max = 2*pi;
else  % if not Grashof, terminate program
  disp('Linkage is not Grashof')
  theta2max = acos((a^2 + d^2 - (b + c)^2)/(2*a*d));
  theta2min = -theta2max;
end

N = 361;   % number of times to perform position calculations
[xB, xC, xP, xQ, xR, xS] = deal(zeros(2, N));          % allocate space for points
[theta2, theta3, theta4] = deal(zeros(1, N));          % allocate space for angles
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
  [ePQ,nPQ] = UnitVector(theta4(i) + gamma2);
  [eRS,nRS] = UnitVector(theta3(i) + gamma4);
 
 % Solve for positions of points B, C and P on the linkage
  xB(:,i) = FindPos(     x0,   a,   e2);
  xC(:,i) = FindPos(     xD,   c,   e4);
  xP(:,i) = FindPos(xC(:,i),   u,   e4);
  xR(:,i) = FindPos(xC(:,i),   v,   e3);
  xQ(:,i) = FindPos(xP(:,i),   y,  ePQ);
  xS(:,i) = FindPos(xR(:,i),   z,  eRS);
end

% Plot paths for points B, C, P, and Q
plot(xB(1,:),xB(2,:), 'Color', [102/255 255/255 255/255])
hold on
plot(xC(1,:),xC(2,:), 'Color', [  0/255 153/255 153/255])
plot(xP(1,:),xP(2,:), 'Color', [204/255 102/255   0/255])
plot(xQ(1,:),xQ(2,:), 'Color', [ 96/255  96/255  96/255])
axis equal
grid on

% Specify angle at which to plot linkage
iTheta = 140;
  
% Plot linkage
plot([x0(1) xB(1,iTheta)],[x0(2) xB(2,iTheta)],'Linewidth',2,'Color','k');
plot([xB(1,iTheta) xC(1,iTheta)],[xB(2,iTheta) xC(2,iTheta)],'Linewidth',2,'Color','k');
plot([xD(1) xC(1,iTheta)],[xD(2) xC(2,iTheta)],'Linewidth',2,'Color','k');
plot([xC(1,iTheta) xP(1,iTheta)],[xC(2,iTheta) xP(2,iTheta)],'Linewidth',2,'Color','k');
plot([xC(1,iTheta) xR(1,iTheta)],[xC(2,iTheta) xR(2,iTheta)],'Linewidth',2,'Color','k');
plot([xP(1,iTheta) xQ(1,iTheta)],[xP(2,iTheta) xQ(2,iTheta)],'Linewidth',2,'Color','k');
plot([xR(1,iTheta) xS(1,iTheta)],[xR(2,iTheta) xS(2,iTheta)],'Linewidth',2,'Color','k');
 
% Plot joints on linkage
plot([x0(1) xD(1) xB(1,iTheta) xC(1,iTheta) xP(1,iTheta) xQ(1,iTheta) xR(1,iTheta) xS(1,iTheta)],...
     [x0(1) xD(2) xB(2,iTheta) xC(2,iTheta) xP(2,iTheta) xQ(2,iTheta) xR(2,iTheta) xS(2,iTheta)],...
     'o','MarkerSize',5,'MarkerFaceColor','k','Color','k');
 
% Plot the labels of each joint
text(       x0(1),       x0(2)-5,'A','HorizontalAlignment','center');
text(       xD(1),       xD(2)-5,'D','HorizontalAlignment','center');
text(xB(1,iTheta),xB(2,iTheta)+5,'B','HorizontalAlignment','center');
text(xC(1,iTheta),xC(2,iTheta)+5,'C','HorizontalAlignment','center');
text(xP(1,iTheta),xP(2,iTheta)+5,'P','HorizontalAlignment','center');
text(xQ(1,iTheta),xQ(2,iTheta)-5,'Q','HorizontalAlignment','center');
text(xR(1,iTheta),xR(2,iTheta)-5,'R','HorizontalAlignment','center');
text(xS(1,iTheta),xS(2,iTheta)+5,'S','HorizontalAlignment','center');
title('Paths of points B, C, P and Q on Claw Mechanism')
xlabel('x-position [mm]')
ylabel('y-position [mm]')
legend('Point B', 'Point C','Point P','Point Q','Location','SouthEast')

% Save the plot
saveas(gcf, 'z_Problem4_20 - plot.png')