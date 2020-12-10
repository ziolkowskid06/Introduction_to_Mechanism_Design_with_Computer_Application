% Conducts a position analysis on the geared fivebar linkage and
% plots the positions of points P and Q
% by Damian Ziolkowski, November 29, 2020

% Prepare Workspace
clear variables; close all; clc;
 
% Linkage dimensions
a = 0.120;           % crank length (m)
b = 0.250;           % coupler 1 length (m)
c = 0.250;           % coupler 2 length (m)
d = 0.180;           % distance between ground pins (m)
u = 0.120;           % length of link on gear 2 (m)
N1 = 24;             % number of teeth on gear 1
N2 = 24;             % number of teeth on gear 2
rho = N1/N2;         % gear ratio
phi = 0;             % offset angle between gears
gamma3 =  20*pi/180; % angle to point P on coupler 1 (CCW rotation)
gamma4 = -20*pi/180; % angle to point Q on coupler 2 (CW rotation)
p = 0.200;           % distance to point P on coupler 1
q = 0.200;           % distance to point Q on coupler 2

% Ground pins
x0 = [0;0];  % ground pin at A (origin)
xD = [d;0];  % ground pin at D
 
N = 361;   % number of times to perform position calculations
[xB,xC,xE_l,xE_r,xP,xQ]       = deal(zeros(2,N));  % allocate space for positions
[theta2,theta3,theta4] = deal(zeros(1,N));  % allocate space for angles
                      % "deal" distributes inputs to outputs

% Perform calculations for each angle
for i = 1:N
  theta2(i) = (i-1)*(2*pi)/(N-1);  % crank angle
  theta5 = -rho*theta2(i) + phi;   % angle of second gear
  [e2,n2] = UnitVector(theta2(i)); % unit vector for crank
  [e5,n5] = UnitVector(theta5);    % unit vector for second gear
  
  xC(:,i) = FindPos(xD, u, e5);    % coords of pin C 
  
  % From now this is fourbar linkage (ACEB)
  dprime = sqrt(xC(1,i)^2 + xC(2,i)^2);  % distance to pin C (ground pins)
  beta = atan2(xC(2,i),xC(1,i));         % angle to pin C (offset)
  r = dprime - a*cos(theta2(i) - beta);
  s = a*sin(theta2(i) - beta);
  f2 = r^2 + s^2;                        % f squared
  delta = acos((b^2+c^2-f2)/(2*b*c));
  g = b - c*cos(delta);
  h = c*sin(delta);
 
  % Calculate rest of fivebar linkage angles 
  theta3(i) = atan2((h*r - g*s),(g*r + h*s)) + beta;
  theta4(i) = theta3(i) + delta;
  
  % Calculate unit vectors
  [e3,n3] = UnitVector(theta3(i));  % unit vector for first coupler
  [e4,n4] = UnitVector(theta4(i));  % unit vector for second coupler
  [eBP,nBP] = UnitVector(theta3(i) + gamma3); % unit vec from B to P
  [eCQ,nCQ] = UnitVector(theta4(i) + gamma4); % unit vec from C to Q
  
  % Solve for positions of points B, E, P and Q on the linkage
  xB(:,i) = FindPos(     x0,   a,   e2);
   % Point E from left and right arm respectively
  xE_l(:,i) = FindPos(xB(:,i),   b,   e3);
  xE_r(:,i) = FindPos(xC(:,i),   c,   e4);
  xP(:,i) = FindPos(xB(:,i),   p,  eBP);
  xQ(:,i) = FindPos(xC(:,i),   q,  eCQ);    
end 

% Plot paths for points E, P and Q
plot(xE_l(1,:),xE_l(2,:), 'Color', [0/255 128/255 255/255])
hold on
plot(xE_r(1,:),xE_r(2,:),'Marker', 'o', 'Color', [155/255 155/255 155/255])
plot(xP(1,:),xP(2,:), 'Color', [204/255 102/255 0/255])
plot(xQ(1,:),xQ(2,:), 'Color', [0/255 153/255 76/255])
axis equal
grid on

% Specify angle at which to plot linkage
iTheta = 120;
 
% Plot the coupler as a triangular patch (filled area)
patch([xB(1,iTheta) xP(1,iTheta) xE_l(1,iTheta)],...
      [xB(2,iTheta) xP(2,iTheta) xE_l(2,iTheta)],[204/255 255/255 229/255]);
patch([xC(1,iTheta) xQ(1,iTheta) xE_r(1,iTheta)],...
      [xC(2,iTheta) xQ(2,iTheta) xE_r(2,iTheta)],[204/255 255/255 229/255]);  
 
% Plot crank and rocker
plot([x0(1) xB(1,iTheta)],[x0(2) xB(2,iTheta)],'Linewidth',2,'Color','k');
plot([xB(1,iTheta) xE_l(1,iTheta)],[xB(2,iTheta) xE_l(2,iTheta)],'Linewidth',2,'Color','k');
plot([xD(1) xC(1,iTheta)],[xD(2) xC(2,iTheta)],'Linewidth',2,'Color','k');
plot([xC(1,iTheta) xE_r(1,iTheta)],[xC(2,iTheta) xE_r(2,iTheta)],'Linewidth',2,'Color','k');
 
% Plot joints on linkage
plot([x0(1) xD(1) xB(1,iTheta) xC(1,iTheta) xE_l(1,iTheta)],...
     [x0(1) xD(2) xB(2,iTheta) xC(2,iTheta) xE_l(2,iTheta)],...
     'o','MarkerSize',5,'MarkerFaceColor','k','Color','k');
 
% Plot the labels of each joint
text(         x0(1),         x0(2)-.015,'A','HorizontalAlignment','center');
text(  xB(1,iTheta),  xB(2,iTheta)+.015,'B','HorizontalAlignment','center');
text(  xC(1,iTheta),  xC(2,iTheta)+.015,'C','HorizontalAlignment','center');
text(         xD(1),         xD(2)-.015,'B','HorizontalAlignment','center');
text(xE_l(1,iTheta),xE_l(2,iTheta)+.015,'E','HorizontalAlignment','center');
text(  xP(1,iTheta),  xP(2,iTheta)+.015,'P','HorizontalAlignment','center');
text(  xQ(1,iTheta),  xQ(2,iTheta)+.015,'Q','HorizontalAlignment','center');
axis equal
grid on
title('Paths of points E, C and P on Fivebar Linkage')
xlabel('x-position [m]')
ylabel('y-position [m]')
legend('Point E - left path','Point E - right path',...
    'Point C', 'Point P','Location','SouthEast')

% Save the plot
saveas(gcf, 'Fivebar_Position_Analysis - plot.png')