% Problem 5.06 (Slider-Crank linkage)
% by Damian Ziolkowski, February 08, 2021

% Prepare Workspace
clear variables; close all; clc;

% Linkage dimensions
a = 0.050;         % crank length (m)
b = 0.200;         % connecting rod length (m)
c = 0.070;          % vertical slider offset (m)
 
% Ground pins
x0 = [0;0];        % ground pin at A (origin)
v0 = [0;0];        % velocity of pin A (zero)

% Angular velocity of crank
omega2 = 10;       % angular velocity of crank (rad/sec)

N = 361;   % number of times to perform position calculations
[xB,xC]           = deal(zeros(2,N));  % allocate space for pins B and C
[vB,vC]           = deal(zeros(2,N));  % allocate space for velocity of B,C
[theta2,theta3,d] = deal(zeros(1,N));  % allocate space for link angles
                 % "deal" distributes inputs to outputs

[omega3,ddot]     = deal(zeros(1,N)); % allocate space for omega3 and ddot

% Perform calculations for each angle
for i = 1:N
  theta2(i) = (i-1)*(2*pi)/(N-1);
  theta3(i) = asin((c - a*sin(theta2(i)))/b);
  d(i) = a*cos(theta2(i)) + b*cos(theta3(i));
 
% Calculate unit vectors
  [e1,n1] = UnitVector(0);
  [e2,n2] = UnitVector(theta2(i));
  [e3,n3] = UnitVector(theta3(i));
  
% Solve for position of point B on the linkage
  xB(:,i) = FindPos(     x0,   a,   e2);
% Not needed but validate the position of point C
  xC(:,i) = FindPos(xB(:,i),   b,   e3);  
  
% Conduct velocity analysis to solve for omega3 and ddot
  A_Mat = [b*n3, -e1];
  b_Vec = -a*omega2*n2;
  omega_Vec = A_Mat\b_Vec;  % solve for velocities
 
  omega3(i) = omega_Vec(1);    % decompose omega_Vec into
  ddot(i)   = omega_Vec(2);    % individual components
  
% Calculate velocity at important points on linkage
  vB(:,i) = FindVel(     v0,   a,    omega2,  n2);    
  vC(:,i) = FindVel(vB(:,i),   b, omega3(i),  n3); 
end

% Code verification - numerical approximation
dt = 2*pi/((N-1)*omega2);             % time increment between calculations
Derivative_Plot(theta2,  theta3,  omega3, dt)   % verify derivatives
Derivative_Plot(theta2, xC(1,:), vC(1,:), dt)   % verify derivatives
 
% Plot the velocity of point C
subplot(1,2,1);
plot(theta2*180/pi,vC(1,:),'Color',[0/255 153/255 76/255])
title('Piston velocity')
xlabel('Crank angle (\circ)')
ylabel('Velocity (m/s)')
grid on
set(gca,'xtick',0:60:360)
xlim([0 360])

subplot(1,2,2);
plot(theta2*180/pi,omega3(:),'Color',[0/255 153/255 76/255])
title('Angular velocity of the connecting rod')
xlabel('Crank angle (\circ)')
ylabel('Angular Velocity (rad/s)')
grid on
set(gca,'xtick',0:60:360)
xlim([0 360])

saveas(gcf, 'z_Problem5_06 - plot.png')
