% Problem 4.6 (ant trajectory)
% by Damian Ziolkowski, December 19, 2020

% Prepare Workspace
clear variables; close all; clc;

% Data
t = 10;                 % time (s)
AB = 200;               % AB length (mm)
BC = 200;               % BC length (mm)
N = 0:0.01:t;           % number of iteration  

% Ground pin
x0 = [0;0];        % ground pin at A (origin)

% Allocate spaces for variables
[xANT,xB,xC] = deal(zeros(2,length(N)));           % allocate space for point B and Ant coordinates
[theta2,theta3, s] = deal(zeros(1,length(N)));     % allocate space for angles and Ant distance

% Main loop
for i = N
   ind = find(N==i); 
   theta2(ind) = 2*i;                  % angle (radians)
   theta3(ind) = -4*i;                 % angle (radians)
   s(ind) = 200-10*i;                  % ant distance (mm)
   
   % Calculate unit vectors
   [eAB, nAB] = UnitVector(theta2(ind));
   [eBC, nBC] = UnitVector(theta3(ind));
   
   % Points and ant coordinates
   xB(:,ind) = FindPos(x0, AB, eAB);
   xC(:,ind) = FindPos(xB(:,ind), BC, eBC);
   xANT(:,ind) = FindPos(xB(:,ind), s(ind), eBC);  
end

% Colors
ant = [  0/255 153/255 76/255];     % ant color
lin = [155/255 155/255 155/255];    % linkage color

% Plot trajectory for the ant
plot(xANT(1,:),xANT(2,:),'Color',ant)
hold on

% Plot bars
plot(    [x0(1) xB(1,end)],    [x0(2) xB(2,end)],'Linewidth',3,'Color',lin);
plot([xB(1,end) xC(1,end)],[xB(2,end) xC(2,end)],'Linewidth',3,'Color',lin);

% Plot joints on linkage
plot([x0(1) xB(1,end) xC(1,end)],...
     [x0(2) xB(2,end) xC(2,end)],...
     'o','MarkerSize',6,'MarkerFaceColor',lin,'Color',lin);

% Plot ant on linkage
plot([xANT(1,end)],...
     [xANT(2,end)],...
     'o','MarkerSize',5,'MarkerFaceColor',ant,'Color',ant);

% Plot the labels of each joint
text(      x0(1)-15,      x0(2)-15,'A','HorizontalAlignment','center');
text(     xB(1,end),  xB(2,end)+30,'B','HorizontalAlignment','center');
text(  xC(1,end)-15,  xC(2,end)-15,'C','HorizontalAlignment','center');
text(xANT(1,end)-10,xANT(2,end)+30,'ANT','HorizontalAlignment','center');
 
 
title('Ant trajectory for 10 seconds')
xlabel('x-position [mm]')
ylabel('y-position [mm]')
axis equal
grid on

% Save the plot
saveas(gcf, 'z_Problem4_06 - plot.png')