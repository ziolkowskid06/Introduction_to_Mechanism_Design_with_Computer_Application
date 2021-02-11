% Problem 6.08 (Threebar Crank-Slider linkage Nonsteady)
% by Damian Ziolkowski, February 11, 2021

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
 
N = 10000;   % number of times to perform position calculations
[xB,xP]           = deal(zeros(2, N)); % allocate space for position of B,P
[vB,vP]           = deal(zeros(2, N)); % allocate space for velocity of B,P
[aB,aP]           = deal(zeros(2, N)); % allocate space for acceleration of B,P
[theta2,theta3]   = deal(zeros(1, N)); % allocate space for link angles
[omega3,alpha3]   = deal(zeros(1, N)); % allocate space for omega3 and alpha3
[b,bdot,bddot]    = deal(zeros(1, N)); % allocate space for b, bdot and bddot
[omega2,alpha2]   = deal(zeros(1, N)); % allocate space for vel and acc of the crank

t = zeros(1,N);      % time (sec)

A = (2*pi*100)/60;    % amp of crank angular accel pulse (rad/sec2)
T = 1;                % period of crank angular acceleration (sec)
B = 2*pi/T;           % freq of crank angular accel pulse (1/sec)
dt = 0.001;            % time increment

% Perform calculation for each angle
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
  b_Vec = -a*omega2(i)*n2;
  omega_Vec = A_Mat\b_Vec;  % solve for velocities
 
  omega3(i) = omega_Vec(1);    % decompose omega_Vec into
  bdot(i)   = omega_Vec(2);    % individual components
 
% Calculate velocity at important points on linkage
  vB(:,i) = FindVel(     v0,   a, omega2(i),  n2);    
  vP(:,i) = FindVel(vB(:,i),   p, omega3(i),  n3);
  
% Conduct acceleration analysis to solve for alpha3 and bddot
  ac = a*omega2(i)^2;        % centripetal acceleration
  at = a*alpha2(i);          % tangential acceleration
  bC = 2*bdot(i)*omega3(i);  % Coriolis acceleration
  bc = b(i)*omega3(i)^2;     % centripetal acceleration
 
  C_Mat = A_Mat;
  d_Vec = -at*n2 + ac*e2 - bC*n3 + bc*e3;
  alpha_Vec = C_Mat\d_Vec;                % solve for angular accelerations
 
  alpha3(i) = alpha_Vec(1);
  bddot(i)  = alpha_Vec(2);  
  
 % Calculate acceleration at important points on linkage
  aB(:,i) = FindAcc(     a0,   a, omega2(i), alpha2(i),   e2,   n2);
  aP(:,i) = FindAcc(aB(:,i),   p, omega3(i), alpha3(i),   e3,   n3);
end

% Plot the acceleration of point P
figure
plot(t,aP(1,:),'Color',[153/255 153/255 153/255])
hold on
plot(t,aP(2,:),'Color',[0/255 153/255 76/255])
legend('a_P_x','a_P_y','location','SouthEast')
title('Acceleration of point P on Threebar Slider-Crank')
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
grid on
grid minor
xlim([0 2])

% Plot crank values
figure
subplot(1,3,1)
plot(t,theta2(1,:),'Color',[0/255 153/255 76/255])
title('Crank Angle')
xlabel('Time (s)')
ylabel('Crank angle (rad)')
grid on
grid minor
xlim([0 2])

subplot(1,3,2)
plot(t,omega2(1,:),'Color',[0/255 153/255 76/255])
title('Angular velocity of crank')
xlabel('Time (s)')
ylabel('Angular velocity (rad/s)')
grid on
grid minor
xlim([0 2])

subplot(1,3,3)
plot(t,alpha2(1,:),'Color',[0/255 153/255 76/255])
title('Angular acceleration pulse of crank')
xlabel('Time (s)')
ylabel('Angular acceleration (m/s^2)')
grid on
grid minor
xlim([0 2])