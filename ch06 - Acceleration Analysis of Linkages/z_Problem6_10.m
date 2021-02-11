% Problem 6.10 (Threebar Slider-Crank linkage)
% by Damian Ziolkowski, February 11, 2021

% Prepare Workspace
clear variables; close all; clc;

% Linkage dimensions
a = 0.030;         % crank length (m)
b = 0.120;         % connecting rod length (m)
c = 0.0;           % vertical slider offset (m)
 
% Ground pins
x0 = [0;0];        % ground pin at A (origin)
v0 = [0;0];        % velocity of pin A (zero)
a0 = [0;0];        % acceleration of pin A (zero)

% Angular velocity of crank
omega2 = (2*pi*1700)/60;     % angular velocity of crank (rad/sec)
alpha2 = 0;                  % angular acceleration of crank (rad/sec2)

N = 361;   % number of times to perform position calculations
[xB,xC]           = deal(zeros(2,N));  % allocate space for pins B and C
[vB,vC]           = deal(zeros(2,N));  % allocate space for velocity of B,C
[aB,aC]           = deal(zeros(2,N));  % allocate space for velocity of B,C
[theta2,theta3] = deal(zeros(1,N));  % allocate space for link angles
                 % "deal" distributes inputs to outputs

[omega3,alpha3]     = deal(zeros(1,N)); % allocate space for omega3 and alpha3
[d,ddot,dddot]      = deal(zeros(1,N)); % allocate space for d, ddot, dddot

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
  
% Conduct acceleration analysis to solve for alpha3 and dddot
  ac = a*omega2^2;           % centripetal acceleration
  at = a*alpha2;             % tangential acceleration
  bc = b*omega3(i)^2;        % centripetal acceleration
 
  C_Mat = A_Mat;
  d_Vec = -at*n2 + ac*e2 + bc*e3;
  alpha_Vec = C_Mat\d_Vec;                % solve for angular accelerations
 
  alpha3(i) = alpha_Vec(1);
  dddot(i)  = alpha_Vec(2);  
  
 % Calculate acceleration at important points on linkage
  aB(:,i) = FindAcc(     a0,   a,    omega2,    alpha2,   e2,   n2);
  aC(:,i) = FindAcc(aB(:,i),   b, omega3(i), alpha3(i),   e3,   n3);
end

% Code verification - numerical approximation
dt = 2*pi/((N-1)*omega2);             % time increment between calculations
Derivative_Plot(theta2,  omega3,  alpha3, dt)   % verify derivatives
Derivative_Plot(theta2,    ddot,   dddot, dt)   % verify derivatives
Derivative_Plot(theta2, vC(1,:), aC(1,:), dt)   % verify derivatives
 
% Plot the acceleration of piston
figure
plot(theta2*180/pi,aC(1,:),'Color',[0/255 153/255 76/255])
title('Acceleration of piston on Slider-Crank')
xlabel('Crank angle (\circ)')
ylabel('Acceleration (m/s^2)')
grid on
set(gca,'xtick',0:60:360)
xlim([0 360])

