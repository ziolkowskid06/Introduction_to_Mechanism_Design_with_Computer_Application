% Problem 4.30 (stamping mechanism)
% by Damian Ziolkowski, December 25, 2020

% Prepare Workspace
clear variables; close all; clc;

% Linkage dimensions
a = 240;                         % AB crank length (mm)
b = 250;                         % BC coupler length (mm)
c = 240;                         % DC rocker length (mm)
d = 220;                         % AD length between ground pins (mm)
u = 210;                         % DP length (mm)
v = 240;                         % handle length (mm)
gamma2 = -50*pi/180;             % CDP angle (rad)
angles = 45:0.1:135;             % handle range (degrees)

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

N = length(angles);   % number of times to perform position calculations
[xB, xC, xP, xZ] = deal(zeros(2, N));                  % allocate space for points
[theta2, theta3, theta4] = deal(zeros(1, N));          % allocate space for angles
                        % "deal" distributes inputs to outputs

% Perform calculation for each angle
for i = 1:N
  theta2(i) = angles(i)*pi/180;
  
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
  [eDP,nDP] = UnitVector(theta4(i) + gamma2);
 
 % Solve for positions of points B, C and P on the linkage
  xB(:,i) = FindPos(     x0,   a,   e2);
  xC(:,i) = FindPos(     xD,   c,   e4);
  xP(:,i) = FindPos(     xD,   u,  eDP);
  xZ(:,i) = FindPos(xB(:,i),   v,   e2);
end

subplot(1,2,1)
% Plot path for point P
plot(xP(1,:),xP(2,:), 'Color', [204/255 102/255   0/255])
hold on
% Specify angle at which to plot linkage
angle = 110;
iTheta = find(angles==110);
  
% Plot linkage
plot([x0(1) xB(1,iTheta)],[x0(2) xB(2,iTheta)],'Linewidth',2,'Color','k');
plot([xB(1,iTheta) xC(1,iTheta)],[xB(2,iTheta) xC(2,iTheta)],'Linewidth',2,'Color','k');
plot([xB(1,iTheta) xZ(1,iTheta)],[xB(2,iTheta) xZ(2,iTheta)],'Linewidth',2,'Color','k');
plot([xD(1) xC(1,iTheta)],[xD(2) xC(2,iTheta)],'Linewidth',2,'Color','k');

patch([xD(1) xC(1,iTheta) xP(1,iTheta)],...
      [xD(2) xC(2,iTheta) xP(2,iTheta)],[200/255 200/255 200/255],...
      'EdgeColor','k','LineWidth',2,'FaceAlpha',0.2); % Make the links semitransparent

% Plot joints on linkage
plot([x0(1) xD(1) xB(1,iTheta) xC(1,iTheta) xP(1,iTheta)],...
     [x0(1) xD(2) xB(2,iTheta) xC(2,iTheta) xP(2,iTheta)],...
     'o','MarkerSize',5,'MarkerFaceColor','k','Color','k');
 
% Plot the labels of each joint
text(x0(1)          ,x0(2)-18      ,'A','HorizontalAlignment','center');
text(xD(1)          ,xD(2)-18      ,'D','HorizontalAlignment','center');
text(xB(1,iTheta)-18,xB(2,iTheta)  ,'B','HorizontalAlignment','center');
text(xC(1,iTheta)   ,xC(2,iTheta)+20,'C','HorizontalAlignment','center');
text(xP(1,iTheta)+18,xP(2,iTheta)+18,'P','HorizontalAlignment','center');
title('Path of point P on Stamping Mechanism')
xlabel('x-position [mm]')
ylabel('y-position [mm]')
grid on

subplot(1,2,2)
plot(angles, xP(2,:), 'Color',[204/255 102/255   0/255])
title('The height of the punch as a function of the handle angle')
xlabel('Handle angle [{\theta}_2]')
ylabel('Height of the punch [mm]')
grid on

% Save the plot
saveas(gcf, 'z_Problem4_30 - plot.png')
