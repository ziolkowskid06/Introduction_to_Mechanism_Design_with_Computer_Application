% Conducts an acceleration analysis on the inverted slider-crank linkage and
% plots the acceleration of point P.
% by Damian Ziolkowski, February 10, 2021
 
% Prepare Workspace
clear variables; close all; clc;
 
% Linkage dimensions
a = 0.080;           % crank length (m)
c = 0.130;           % rocker length (m)
d = 0.200;           % length between ground pins (m)
p = 0.350;           % slider length (m)
delta = 60*pi/180;   % angle between slider and rocker (converted to rad)
eps = 0.000001;      % tiny number added to theta2 to keep it in bounds
h = c*sin(delta);    % h is a constant, only calculate it once
 
% Ground pins
x0 = [0;0];   % ground pin at A (origin)
xD = [d;0];   % ground pin at D
v0 = [0;0];    % velocity of pin A (zero)
a0 = [0;0];    % acceleration of pin A (zero)

% Grashof check
if ((a + d) < c*sin(delta))
  disp('Linkage cannot be assembled')
  return
else
  if (abs(a - d) >= c*sin(delta))
    disp('Linkage is Grashof')
    theta2min = 0;
    theta2max = 2*pi;
  else
    disp('Linkage is not Grashof')
    theta2min = acos((a^2 + d^2 - (c*sin(delta))^2)/(2*a*d)) + eps;
    theta2max = 2*pi - theta2min;
  end
end

% Angular velocity and acceleration of crank
omega2 = 10;      % angular velocity of crank (rad/s)
alpha2 = 0;       % angular acceleration of crank (rad/s^2)

N = 361;   % number of times to perform position calculations
[xB, xC, xP] = deal(zeros(2, N));          % allocate space for positions
[vB, vC, vP] = deal(zeros(2, N));          % allocate space for velocities
[aB, aC, aP] = deal(zeros(2, N));          % allocate space for accelerations
[theta2, theta3, theta4] = deal(zeros(1, N)); % allocate space for angles
                        % "deal" distributes inputs to outputs
[omega3,alpha3] = deal(zeros(1,N));   % allocate space for omega3 and alpha3
[b,bdot,bddot]  = deal(zeros(1,N));   % allocate space for b, bdot, bddot

% Perform calculations for each angle
for i = 1:N
  theta2(i) = (i-1)*(theta2max - theta2min)/(N-1) + theta2min;
  r = d - a*cos(theta2(i));
  s = a*sin(theta2(i));
  f2 = r^2 + s^2;  % f squared
 
  b(i) = c * cos(delta) + sqrt(f2 - h^2);
  g = b(i) - c*cos(delta);
 
  theta3(i) = atan2((h*r - g*s),(g*r + h*s));
  theta4(i) = theta3(i) + delta;
 
 % Calculate unit vectors
  [e2,n2]   = UnitVector(theta2(i));
  [e3,n3]   = UnitVector(theta3(i));
  [e4,n4]   = UnitVector(theta4(i));
 
 % Solve for positions of points B, C and P on the linkage
  xB(:,i) = FindPos(     x0,a,e2);
  xC(:,i) = FindPos(     xD,c,e4);
  xP(:,i) = FindPos(xB(:,i),p,e3);
  
  % Conduct velocity analysis to solve for omega3 and bdot
  A_Mat = [b(i)*n3-c*n4, e3];
  b_Vec = -a*omega2*n2;
  omega_Vec = A_Mat\b_Vec;  % solve for angular velocities
 
  omega3(i) = omega_Vec(1);    % decompose omega_Vec into
    bdot(i) = omega_Vec(2);    % individual components
     omega4 = omega3; 
 
% Calculate velocity at important points on linkage
  vB(:,i) = FindVel(     v0,   a,    omega2,   n2);
  vC(:,i) = FindVel(     v0,   c, omega4(i),   n4);
  vP(:,i) = FindVel(vB(:,i),   p, omega3(i),   n3); 

  % Conduct acceleration analysis to solve for alpha3 and bddot
  ac = a*omega2^2;             % centripetal acceleration
  at = a*alpha2;               % tangential acceleration
  bc = b(i)*omega3(i)^2;       % centripetal acceleration
  bC = 2*bdot(i)*omega3(i);    % Coriolis acceleration
  cc = c*omega3(i)^2;          % centripetal acceleration
 
  C_Mat = A_Mat;
  d_Vec = -at*n2 + ac*e2 - bC*n3 + bc*e3 - cc*e4;
  alpha_Vec = C_Mat\d_Vec;  % solve for angular accelerations
 
  alpha3(i) = alpha_Vec(1);
  bddot(i)  = alpha_Vec(2);
 
% Find acceleration of pins
  aB(:,i) = FindAcc(     a0, a,    omega2,    alpha2,  e2,  n2);
  aC(:,i) = FindAcc(     a0, c, omega4(i), alpha3(i),  e4,  n4);
  aP(:,i) = FindAcc(aB(:,i), p, omega3(i), alpha3(i),  e3,  n3);
end

% Code verification - numerical approximation
dt = 2*pi/((N-1)*omega2);             % time increment between calculations
Derivative_Plot(theta2, omega3, alpha3, dt)   % verify derivatives
Derivative_Plot(theta2,   bdot,  bddot, dt)   % verify derivatives
Derivative_Plot(theta2,vP(1,:),aP(1,:), dt)   % verify derivatives
 
% Plot the acceleration of point P
figure
plot(theta2*180/pi,aP(1,:),'Color',[153/255 153/255 153/255])
hold on
plot(theta2*180/pi,aP(2,:),'Color',[0/255 153/255 76/255])
legend('a_P_x','a_P_y','Location','Southeast')
title('Acceleration of point P on Inverted Slider-Crank Linkage')
xlabel('Crank angle (\circ)')
ylabel('Acceleration (m/s^2)')
grid on
set(gca,'xtick',0:60:360)
xlim([0 360])

% Save the plot
saveas(gcf, 'InvSlider_Acceleration_Analysis - plot.png')