% Conducts a position analysis of the Stephenson Type I sixbar linkage
% and plots the positions of points E, F and G.
% by Damian Ziolkowski, December 11, 2020
 
% Prepare Workspace
clear variables; close all; clc;
 
% Linkage dimensions
a = 0.070;           % crank length (m)
b = 0.100;           % coupler length (m)
c = 0.090;           % rocker length (m)
d = 0.110;           % length between ground pins (m)
p = 0.150;           % length to third pin on crank triangle (m)
q = 0.150;           % length to third pin on rocker triangle (m)
u = 0.120;           % length of link 5 (m)
v = 0.160;           % length of link 6 (m)
gamma2 = 20*pi/180;  % internal angle of crank triangle (CCW rotation)
gamma4 = -20*pi/180; % internal angle of rocker triangle (CW rotation)

% Ground pins
x0 = [ 0; 0];    % ground pin at A (origin)
xD = [ d; 0];    % ground pin at D
 
% Allocate space for variables
N = 361;   % number of times to perform position calculations
[xB,xC,xE,xF,xG]                     = deal(zeros(2,N)); % points
[theta2,theta3,theta4,theta5,theta6] = deal(zeros(1,N)); % angles
                                    % "deal" distributes inputs to outputs

% Perform calculations for every angle
for i = 1:N
 
  % Solve lower fourbar linkage
  theta2(i) = (i-1)*(2*pi)/(N-1);      % crank angle
  r = d - a*cos(theta2(i));
  s = a*sin(theta2(i));
  f2 = r^2 + s^2;                      % f squared
  delta = acos((b^2+c^2-f2)/(2*b*c));  % angle between coupler and rocker
  g = b - c*cos(delta);
  h = c*sin(delta);
  theta3(i) = atan2((h*r - g*s),(g*r + h*s)); % coupler angle
  theta4(i) = theta3(i) + delta;              % rocker angle
  
  % Calculate unit vectors
  [e2,n2] = UnitVector(theta2(i));
  [e3,n3] = UnitVector(theta3(i));
  [e4,n4] = UnitVector(theta4(i));
  [eAE,nAE] = UnitVector(theta2(i) + gamma2);
  [eDF,nDF] = UnitVector(theta4(i) + gamma4);
 
  % Solve for positions of points B, C, E, F
  xB(:,i) = FindPos(x0, a,  e2);
  xC(:,i) = FindPos(xD, c,  e4);
  xE(:,i) = FindPos(x0, p, eAE);
  xF(:,i) = FindPos(xD, q, eDF);
  
  % Solve upper fourbar linkage
  xFB = xF(1,i) - xB(1,i);   yFB = xF(2,i) - xB(2,i);  
  xEB = xE(1,i) - xB(1,i);   yEB = xE(2,i) - xB(2,i);
  beta =  atan2(yFB, xFB);
  alpha = atan2(yEB, xEB);
  aPrime = sqrt(xEB^2 + yEB^2);        % virtual crank length on upper fourbar
  dPrime = sqrt(xFB^2 + yFB^2);        % virtual ground length on upper fourbar
  theta2Prime = alpha - beta;          % virtual crank angle on upper fourbar
  r = dPrime - aPrime*cos(theta2Prime);
  s = aPrime*sin(theta2Prime);
  f2 = r^2 + s^2;
  delta = acos((u^2+v^2-f2)/(2*u*v));
  g = u - v*cos(delta);
  h = v*sin(delta);
  theta5Prime = atan2((h*r - g*s),(g*r + h*s));   % coupler and rocker 
  theta6Prime = theta5Prime + delta;              % angles on upper fourbar
  
  % Return angles to fixed coordinate system
  theta5(i) = theta5Prime + beta;                 
  theta6(i) = theta6Prime + beta;                
 
  % Calculate remaining unit vectors
  [e5,n5] = UnitVector(theta5(i));
  [e6,n6] = UnitVector(theta6(i));
                 
  % Calculate position of point G
  xG(:,i) = FindPos(xE(:,i), u, e5); 
end

% Give shades to plot
cGre = DefineColor([  0 148  74]); % greenscale
cBlk = DefineColor([  0   0   0]); % grayscale

% Plot trajectories of points E, F, and G 
figure; hold on
plot(xE(1,:), xE(2,:), 'Color', cGre(1,:))
plot(xF(1,:), xF(2,:), 'Color', cGre(7,:))
plot(xG(1,:), xG(2,:), 'Color', cBlk(5,:))

% Specify angle at which to plot linkage
iTheta = 120;  

% Plot the two three-pin links as triangular patches
patch([x0(1) xB(1,iTheta) xE(1,iTheta)],...
      [x0(2) xB(2,iTheta) xE(2,iTheta)],cGre(9,:),...
      'EdgeColor',cBlk(6,:),'LineWidth',2,'FaceAlpha',0.5); % Make the links semitransparent
patch([xD(1) xF(1,iTheta) xC(1,iTheta)],...
      [xD(2) xF(2,iTheta) xC(2,iTheta)],cGre(10,:),...
      'EdgeColor',cBlk(6,:),'LineWidth',2,'FaceAlpha',0.5); % Make the links semitransparent

% Plot the two-pin links
plot([xB(1,iTheta) xC(1,iTheta)],[xB(2,iTheta) xC(2,iTheta)],...
     'Linewidth',2,'Color',cBlk(5,:));
plot([xE(1,iTheta) xG(1,iTheta)],[xE(2,iTheta) xG(2,iTheta)],...
     'Linewidth',2,'Color',cBlk(5,:));
plot([xF(1,iTheta) xG(1,iTheta)],[xF(2,iTheta) xG(2,iTheta)],...
     'Linewidth',2,'Color',cBlk(5,:));
 
% Plot fixed pins
 plot([x0(1) xD(1)],[x0(2) xD(2)],'o','MarkerSize',5,...
      'MarkerFaceColor',cBlk(1,:),'Color',cBlk(1,:));
 
% Plot joints on linkage
 plot([xB(1,iTheta) xC(1,iTheta) xE(1,iTheta) xF(1,iTheta) xG(1,iTheta)],...
      [xB(2,iTheta) xC(2,iTheta) xE(2,iTheta) xF(2,iTheta) xG(2,iTheta)],...
      'o','MarkerSize',5,'MarkerFaceColor',cBlk(6,:),'Color',cBlk(6,:));

% Plot the labels of each joint
text(x0(1),x0(2)-0.012,'A','HorizontalAlignment','center');
text(xD(1),xD(2)-0.012,'D','HorizontalAlignment','center');
text(xB(1,iTheta),xB(2,iTheta)+0.012,'B','HorizontalAlignment','center');
text(xC(1,iTheta)+0.012,xC(2,iTheta),'C','HorizontalAlignment','center');
text(xE(1,iTheta)-0.012,xE(2,iTheta),'E','HorizontalAlignment','center');
text(xF(1,iTheta),xF(2,iTheta)+0.012,'F','HorizontalAlignment','center');
text(xG(1,iTheta),xG(2,iTheta)+0.012,'G','HorizontalAlignment','center');

% Plot parameters
axis equal
grid on
title('Paths of points E, F and G on Stephenson Type I sixbar linkage')
xlabel('x-position [m]')
ylabel('y-position [m]')
legend('Point E', 'Point F', 'Point G','Location','SouthEast')

% Save the plot
saveas(gcf, 'Sixbar_StephensonType1_Position_Analysis - plot.png')
