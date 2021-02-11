% Problem 6.07 (Threebar Crank-Slider linkage)
% by Damian Ziolkowski, February 10, 2021

% Prepare Workspace
clear variables; close all; clc;
 
% Linkage dimensions
a = 2;         % crank length (m)
c = 3;         % y-distance between ground pins (m)
d = 4;         % x-distance between ground pins (m)
p = 5;         % slider length (m)
 
% Ground pins
xA = [0;0];        % point A (the origin)
xD = [d;c];        % point D
v0 = [0;0];        % velocity of pin A (zero)
a0 = [0;0];        % acceleration of pin A (zero)

% Angular velocity of crank
n = 100;                % rotations of crank  (rpm) 
omega2 = 2*pi/60;       % angular velocity of crank (rad/sec)
alpha2 = 0;       % angular acceleration of crank (rad/sec^2)

 
N = 361;   % number of times to perform position calculations
[xB,xP]           = deal(zeros(2, N)); % allocate space for position of B,P
[vB,vP]           = deal(zeros(2, N)); % allocate space for velocity of B,P
[aB,aP]           = deal(zeros(2, N)); % allocate space for acceleration of B,P
[theta2,theta3]   = deal(zeros(1, N)); % allocate space for link angles
[omega3,alpha3]   = deal(zeros(1, N)); % allocate space for omega3 and alpha3
[b,bdot,bddot]    = deal(zeros(1, N)); % allocate space for b, bdot and bddot
 
% Perform calculation for each angle
for i = 1:N
  % Map iterations to angles to make sure crank makes complete rotation  
  theta2(i) = (i-1)*(2*pi)/(N-1);
  % Use "atan2" to calculate angle for all four quadrants
  theta3(i) = atan2(-a*sin(theta2(i)) + c, d - a*cos(theta2(i)));
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
  
% Conduct acceleration analysis to solve for alpha3 and bddot
  ac = a*omega2^2;           % centripetal acceleration
  at = a*alpha2;             % tangential acceleration
  bC = 2*bdot(i)*omega3(i);  % Coriolis acceleration
  bc = b(i)*omega3(i)^2;     % centripetal acceleration
 
  C_Mat = A_Mat;
  d_Vec = -at*n2 + ac*e2 - bC*n3 + bc*e3;
  alpha_Vec = C_Mat\d_Vec;                % solve for angular accelerations
 
  alpha3(i) = alpha_Vec(1);
  bddot(i)  = alpha_Vec(2);  
  
 % Calculate acceleration at important points on linkage
  aB(:,i) = FindAcc(     a0,   a,    omega2,    alpha2,   e2,   n2);
  aP(:,i) = FindAcc(aB(:,i),   p, omega3(i), alpha3(i),   e3,   n3);
end

% Code verification - numerical approximation
dt = 2*pi/((N-1)*omega2);             % time increment between calculations
Derivative_Plot(theta2,  omega3,  alpha3, dt)   % verify derivatives
Derivative_Plot(theta2,    bdot,   bddot, dt)   % verify derivatives
Derivative_Plot(theta2, vP(1,:), aP(1,:), dt)   % verify derivatives
 
% Plot the acceleration of point P
figure
plot(theta2*180/pi,aP(1,:),'Color',[153/255 153/255 153/255])
hold on
plot(theta2*180/pi,aP(2,:),'Color',[0/255 153/255 76/255])
legend('a_P_x','a_P_y','location','SouthEast')
title('Acceleration of point P on Threebar Slider-Crank')
xlabel('Crank angle (\circ)')
ylabel('Acceleration (m/s^2)')
grid on
set(gca,'xtick',0:60:360)
xlim([0 360])

% Save the plot
saveas(gcf, 'z_Problem6_07 - plot.png')