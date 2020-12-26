% Conducts a position analysis of the Stephenson Type III sixbar linkage
% and plots the positions of points C, E and G.
% by Damian Ziolkowski, December 12, 2020

% Prepare Workspace
clear variables; close all; clc;
 
% Linkage dimensions
a = 0.080;           % AB crank length (m)
b = 0.105;           % BC coupler length (m)
c = 0.120;           % CD rocker length (m)
d = 0.135;           % AD length between ground pins (m)
p = 0.225;           % BE length on coupler (m)
q = 0.195;           % AF length on ground (m)
u = 0.180;           % EG length of link 5 (m)
v = 0.180;           % FG length of link 6 (m)
gamma1 = -20*pi/180; % internal angle of ground (CW rotation)
gamma3 = 20*pi/180;  % internal angle of coupler triangle (CCW rotation)

% Ground pins
x0 = [ 0; 0];    % ground pin at A (origin)
xD = [ d; 0];    % ground pin at D
[eAF,nAF] = UnitVector(gamma1);  
xF = FindPos(x0, q, eAF);   % ground pin at F

% Allocate space for variables
N = 361;   % number of times to perform position calculations
[xB,xC,xE,xG]                        = deal(zeros(2,N)); % points
[theta2,theta3,theta4,theta5,theta6] = deal(zeros(1,N)); % angles
                                    % "deal" distributes inputs to outputs

% Perform calculations for every angle
for i = 1:N
 
  % Solve left fourbar linkage
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
  [eBE,nBE] = UnitVector(theta3(i) + gamma3);
  
  % Solve for positions of points B, C, E
  xB(:,i) = FindPos(x0, a, e2);
  xC(:,i) = FindPos(xB(:,i), b, e3);
  xE(:,i) = FindPos(xB(:,i), p, eBE);
  
  % Solve right fourbar linkage
  xFC = xF(1) - xC(1,i);   
  yFC = xF(2) - xC(2,i);  
  xEC = xE(1,i) - xC(1,i);   
  yEC = xE(2,i) - xC(2,i);
  beta =  atan2(yFC, xFC);
  alpha = atan2(yEC, xEC);
  theta2Prime = alpha - beta;          % virtual crank angle on upper fourbar
  aPrime = sqrt(xEC^2 + yEC^2);        % virtual crank length on right fourbar
  dPrime = sqrt(xFC^2 + yFC^2);        % virtual ground length on upper fourbar
 
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

% Plot trajectories of points C, E, and G 
figure; hold on
plot(xE(1,:), xE(2,:), 'Color', cGre(1,:))
plot(xC(1,:), xC(2,:), 'Color', cGre(7,:))
plot(xG(1,:), xG(2,:), 'Color', cBlk(5,:))

% Specify angle at which to plot linkage
iTheta = 120;  

% Plot the two three-pin links as triangular patches
patch([xB(1,iTheta) xC(1,iTheta) xE(1,iTheta)],...
      [xB(2,iTheta) xC(2,iTheta) xE(2,iTheta)],cGre(9,:),...
      'EdgeColor',cBlk(6,:),'LineWidth',2,'FaceAlpha',0.5); % Make the links semitransparent
patch([x0(1) xD(1) xF(1)],...
      [x0(2) xD(2) xF(2)],cBlk(10,:),...
      'EdgeColor',cBlk(6,:),'LineWidth',2,'FaceAlpha',0.5); % Make the links semitransparent

% Plot the two-pin links
plot([x0(1) xB(1,iTheta)],[x0(2) xB(2,iTheta)],...
     'Linewidth',2,'Color',cBlk(5,:));
plot([xC(1,iTheta) xD(1)],[xC(2,iTheta) xD(2)],...
     'Linewidth',2,'Color',cBlk(5,:));
plot([xF(1) xG(1,iTheta)],[xF(2) xG(2,iTheta)],...
     'Linewidth',2,'Color',cBlk(5,:));
plot([xE(1,iTheta) xG(1,iTheta)],[xE(2,iTheta) xG(2,iTheta)],...
     'Linewidth',2,'Color',cBlk(5,:));
 
% Plot fixed pins
 plot([x0(1) xD(1) xF(1)],[x0(2) xD(2) xF(2)],'o','MarkerSize',5,...
      'MarkerFaceColor',cBlk(1,:),'Color',cBlk(1,:));
 
% Plot joints on linkage
 plot([xB(1,iTheta) xC(1,iTheta) xE(1,iTheta) xG(1,iTheta)],...
      [xB(2,iTheta) xC(2,iTheta) xE(2,iTheta) xG(2,iTheta)],...
      'o','MarkerSize',5,'MarkerFaceColor',cBlk(6,:),'Color',cBlk(6,:));

% Plot the labels of each joint
text(x0(1)+0.005,x0(2)+0.012,'A','HorizontalAlignment','center');
text(xD(1)+0.005,xD(2)+0.012,'D','HorizontalAlignment','center');
text(xB(1,iTheta),xB(2,iTheta)+0.012,'B','HorizontalAlignment','center');
text(xC(1,iTheta)+0.012,xC(2,iTheta),'C','HorizontalAlignment','center');
text(xE(1,iTheta)-0.005,xE(2,iTheta)+0.012,'E','HorizontalAlignment','center');
text(xF(1),xF(2)+0.012,'F','HorizontalAlignment','center');
text(xG(1,iTheta),xG(2,iTheta)+0.012,'G','HorizontalAlignment','center');

% Plot parameters
axis equal
grid on
title('Paths of points C, E and G on Stephenson Type III sixbar linkage')
xlabel('x-position [m]')
ylabel('y-position [m]')
legend('Point E', 'Point C', 'Point G','Location','SouthWest')

% Save the plot
saveas(gcf, 'Sixbar_StephensonType3_Position_Analysis - plot.png')
