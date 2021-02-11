% Problem 6.11 (Threebar Slider-Crank linkage Nonsteady)
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

N = 1000;   % number of time steps to calculate
[xB,xC]           = deal(zeros(2,N));  % allocate space for pins B and C
[vB,vC]           = deal(zeros(2,N));  % allocate space for velocity of B,C
[aB,aC]           = deal(zeros(2,N));  % allocate space for velocity of B,C
[theta2,theta3] = deal(zeros(1,N));  % allocate space for link angles
                 % "deal" distributes inputs to outputs

[omega3,alpha3] = deal(zeros(1,N)); % allocate space for omega3 and alpha3
[d,ddot,dddot]  = deal(zeros(1,N)); % allocate space for d, ddot, dddot
[omega2,alpha2] = deal(zeros(1,N)); % angular velocity and acceleration of the crank
t = zeros(1,N);      % time (sec)

A = (2*pi*500)/(2*60);   % amp of crank angular accel pulse (rad/sec2)
T = 2;                   % period of crank angular acceleration (sec)
B = 2*pi/T;              % freq of crank angular accel pulse (1/sec)
dt = 0.01;               % time increment

% Perform calculations for each angle
for i = 1:N
  t(i) = i*dt;       % calculate time
  % Determine current values of crank angular acceleration, velocity and position
  if (t(i) <= T)  % calculate crank angle, vel, accel for t < T
    alpha2(i) = A*(1-cos(2*pi*t(i)/T));
    omega2(i) = A*t(i) - (1/B)*sin(B*t(i));
    theta2(i) = A*t(i)^2/2 + (1/(B^2))*(cos(B*t(i))-1);
  else            % calculate crank angle, vel, accel for t > T
    alpha2(i) = 0;
    omega2(i) = A*T;
    theta2(i) = A*T*t(i) - T^2*A/2;
  end
  
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
  b_Vec = -a*omega2(i)*n2;
  omega_Vec = A_Mat\b_Vec;  % solve for velocities
 
  omega3(i) = omega_Vec(1);    % decompose omega_Vec into
  ddot(i)   = omega_Vec(2);    % individual components
  
% Calculate velocity at important points on linkage
  vB(:,i) = FindVel(     v0,   a, omega2(i),  n2);    
  vC(:,i) = FindVel(vB(:,i),   b, omega3(i),  n3); 
  
% Conduct acceleration analysis to solve for alpha3 and dddot
  ac = a*omega2(i)^2;           % centripetal acceleration
  at = a*alpha2(i);             % tangential acceleration
  bc = b*omega3(i)^2;        % centripetal acceleration
 
  C_Mat = A_Mat;
  d_Vec = -at*n2 + ac*e2 + bc*e3;
  alpha_Vec = C_Mat\d_Vec;                % solve for angular accelerations
 
  alpha3(i) = alpha_Vec(1);
  dddot(i)  = alpha_Vec(2);  
  
 % Calculate acceleration at important points on linkage
  aB(:,i) = FindAcc(     a0,   a, omega2(i), alpha2(i),   e2,   n2);
  aC(:,i) = FindAcc(aB(:,i),   b, omega3(i), alpha3(i),   e3,   n3);
end
 
% Plot piston values
figure
subplot(1,3,1)
plot(t,xC(1,:),'Color',[0/255 153/255 76/255])
title('Position of piston on slider-crank linkage')
xlabel('Time (s)')
ylabel('Position (m)')
grid on
grid minor
xlim([0 3])

subplot(1,3,2)
plot(t,vC(1,:),'Color',[0/255 153/255 76/255])
title('Velocity of piston on slider-crank linkage')
xlabel('Time (s)')
ylabel('Velocity (m/s)')
grid on
grid minor
xlim([0 3])

subplot(1,3,3)
plot(t,aC(1,:),'Color',[0/255 153/255 76/255])
title('Acceleration of piston on slider-crank linkage')
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
grid on
grid minor
xlim([0 3])

% Plot crank values
figure
subplot(1,3,1)
plot(t,theta2(1,:),'Color',[0/255 153/255 76/255])
title('Crank Angle')
xlabel('Time (s)')
ylabel('Crank angle (rad)')
grid on
grid minor
xlim([0 3])

subplot(1,3,2)
plot(t,omega2(1,:),'Color',[0/255 153/255 76/255])
title('Angular velocity of crank')
xlabel('Time (s)')
ylabel('Angular velocity (rad/s)')
grid on
grid minor
xlim([0 3])

subplot(1,3,3)
plot(t,alpha2(1,:),'Color',[0/255 153/255 76/255])
title('Angular acceleration pulse of crank')
xlabel('Time (s)')
ylabel('Angular acceleration (m/s^2)')
grid on
grid minor
xlim([0 3])