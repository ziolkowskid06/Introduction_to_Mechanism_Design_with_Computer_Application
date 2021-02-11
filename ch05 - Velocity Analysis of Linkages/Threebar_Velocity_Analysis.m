% Conducts a velocity analysis on the threebar crank-slider linkage
% and plots the velocity of the point P
% by Damian Ziolkowski, February 5, 2021

% Prepare Workspace
clear variables; close all; clc;
 
% Linkage dimensions
a = 0.100;         % crank length (m)
d = 0.150;         % length between ground pins (m)
p = 0.300;         % slider length (m)
 
% Ground pins
xA = [0;0];        % point A (the origin)
xD = [d;0];        % point D
v0 = [0;0];        % velocity of pin A (zero)

% Angular velocity of crank
omega2 = 10;       % angular velocity of crank (rad/sec)
 
N = 361;   % number of times to perform position calculations
[xB,xP]           = deal(zeros(2, N)); % allocate space for position of B,P
[vB,vP]           = deal(zeros(2, N)); % allocate space for velocity of B,P
[theta2,theta3,b] = deal(zeros(1, N)); % allocate space for link angles
[omega3,bdot]     = deal(zeros(1, N)); % allocate space for omega3 and bdot
 
% Perform calculation for each angle
for i = 1:N
  % Map iterations to angles to make sure crank makes complete rotation  
  theta2(i) = (i-1)*(2*pi)/(N-1);
  % Use "atan2" to calculate angle for all four quadrants
  theta3(i) = atan2(-a*sin(theta2(i)), d - a*cos(theta2(i)));
  b(i) = (d - a*cos(theta2(i)))/cos(theta3(i));
 
% Calculate unit vectors
  [e1,n1]   = UnitVector(0);
  [e2,n2]   = UnitVector(theta2(i));
  [e3,n3]   = UnitVector(theta3(i));
 
% Solve for positions of points B and P on the linkage
  xB(:,i) = FindPos(     xA, a, e2);
  xP(:,i) = FindPos(xB(:,i), p, e3);
  
% Conduct velocity analysis to solve for omega3 and bdot
  A_Mat = [b(i)*n3 e3];
  b_Vec = -a*omega2*n2;
  omega_Vec = A_Mat\b_Vec;  % solve for velocities
 
  omega3(i) = omega_Vec(1);    % decompose omega_Vec into
  bdot(i)   = omega_Vec(2);    % individual components
 
% Calculate velocity at important points on linkage
  vB(:,i) = FindVel(     v0,   a,    omega2,  n2);    
  vP(:,i) = FindVel(vB(:,i),   p, omega3(i),  n3);
end

% Code verification - numerical approximation
dt = 2*pi/((N-1)*omega2);             % time increment between calculations
Derivative_Plot(theta2, theta3, omega3, dt)   % verify derivatives
Derivative_Plot(theta2, b, bdot, dt)          % verify derivatives
Derivative_Plot(theta2, xP(1,:), vP(1,:), dt) % verify derivatives
 
% Plot the velocity of point P
figure
plot(theta2*180/pi,vP(1,:),'Color',[153/255 153/255 153/255])
hold on
plot(theta2*180/pi,vP(2,:),'Color',[0/255 153/255 76/255])
legend('v_P_x','v_P_y','location','SouthEast')
title('Velocity of point P on Threebar Slider-Crank')
xlabel('Crank angle (\circ)')
ylabel('Velocity (m/s)')
grid on
set(gca,'xtick',0:60:360)
xlim([0 360])

% Save the plot
saveas(gcf, 'Threebar_Velocity_Analysis - plot1.png')