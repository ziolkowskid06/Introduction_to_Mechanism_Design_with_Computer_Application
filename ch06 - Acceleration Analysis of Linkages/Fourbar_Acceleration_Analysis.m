% Conducts an acceleration analysis on the fourbar linkage and plots the
% acceleration of point P.
% by Damian Ziolkowski, February 10, 2021
 
% Prepare Workspace
clear variables; close all; clc;

% Linkage dimensions
a = 0.130;         % crank length (m)
b = 0.200;         % coupler length (m)
c = 0.170;         % rocker length (m)
d = 0.220;         % length between ground pins (m)
p = 0.150;         % length from B to P (m)
gamma = 20*pi/180; % angle between BP and coupler (converted to rad)
 
% Ground pins
x0 = [0; 0];   % ground pin at A (origin)
xD = [d; 0];   % ground pin at D
v0 = [0;0];    % velocity of pin A (zero)
a0 = [0;0];    % acceleration of pin A (zero)

% Angular velocity and acceleration of crank
omega2 = 10;      % angular velocity of crank (rad/s)
alpha2 = 0;       % angular acceleration of crank (rad/s^2)

% Grashof Check
S = min([a b c d]);  % length of shortest link
L = max([a b c d]);  % length of longest link
T = sum([a b c d]);  % total of all link lengths
PQ = T - S - L;  % length of P plus length of Q
if (S+L < PQ)  % Grashof condition
  disp('Linkage is Grashof.')
  theta2min = 0;
  theta2max = 2*pi;
else  % if not Grashof, terminate program
  disp('Linkage is not Grashof')
  theta2max = acos((a^2 + d^2 - (b + c)^2)/(2*a*d));
  theta2min = -theta2max;
end

N = 361;   % number of times to perform position calculations
[xB, xC, xP] = deal(zeros(2, N));      % allocate space for points
[vB, vC, vP] = deal(zeros(2,N));       % allocate space for vel of B, C, P
[aB, aC, aP] = deal(zeros(2,N));       % allocate space for acc of B, C, P
[theta2, theta3, theta4] = deal(zeros(1, N)); % allocate space for angles
[omega3,omega4] = deal(zeros(1,N));  % allocate for angular velocity
[alpha3,alpha4] = deal(zeros(1,N));  % allocate for angular acceleration

% Perform calculations for each angle
for i = 1:N
  theta2(i) = (i-1)*(theta2max - theta2min)/(N-1) + theta2min;  % crank angle
 
 % Conduct position analysis to solve for theta3 and theta4
  r = d - a*cos(theta2(i));
  s = a*sin(theta2(i));
  f2 = r^2 + s^2;                     % f squared
  delta = acos((b^2+c^2-f2)/(2*b*c)); % angle between coupler and rocker
 
  g = b - c*cos(delta);
  h = c*sin(delta);
 
  theta3(i) = atan2((h*r - g*s),(g*r + h*s));  % coupler angle
  theta4(i) = theta3(i) + delta;               % rocker angle
 
 % Calculate unit vectors
  [e1,n1] = UnitVector(0);
  [e2,n2] = UnitVector(theta2(i));
  [e3,n3] = UnitVector(theta3(i));
  [e4,n4] = UnitVector(theta4(i));
  [eBP,nBP] = UnitVector(theta3(i) + gamma);
 
 % Solve for positions of points B, C and P on the linkage
  xB(:,i) = FindPos(     x0,   a,   e2);
  xC(:,i) = FindPos(     xD,   c,   e4);
  xP(:,i) = FindPos(xB(:,i),   p,  eBP);
  
 % Conduct velocity analysis to solve for omega3 and omega4
  A_Mat = [b*n3    -c*n4];
  b_Vec = -a*omega2*n2;
  omega_Vec = A_Mat\b_Vec;  % solve for angular velocities
 
  omega3(i) = omega_Vec(1);    % decompose omega_Vec into
  omega4(i) = omega_Vec(2);    % individual components
 
% Calculate velocity at important points on linkage
  vB(:,i) = FindVel(     v0,   a,    omega2,   n2);
  vC(:,i) = FindVel(     v0,   c, omega4(i),   n4);
  vP(:,i) = FindVel(vB(:,i),   p, omega3(i),  nBP); 
  
% Conduct acceleration analysis to solve for alpha3 and alpha4
  ac = a*omega2^2;       % centripetal accel |AB|
  bc = b*omega3(i)^2;    % centripetal accel |BC|
  cc = c*omega4(i)^2;    % centripetal accel |CD|
  pc = p*omega3(i)^2;    % centripetal accel |BP|
 
  C_Mat = A_Mat;
  d_Vec = ac*e2 + bc*e3 - cc*e4;
  alpha_Vec = C_Mat\d_Vec;  % solve for angular accelerations
 
  alpha3(i) = alpha_Vec(1);
  alpha4(i) = alpha_Vec(2);
 
% Find acceleration of pins
  aB(:,i) = FindAcc(     a0, a,    omega2,    alpha2,  e2,  n2);
  aC(:,i) = FindAcc(     a0, c, omega4(i), alpha4(i),  e4,  n4);
  aP(:,i) = FindAcc(aB(:,i), p, omega3(i), alpha3(i), eBP, nBP);
end

% Code verification - numerical approximation
dt = 2*pi/((N-1)*omega2);             % time increment between calculations
Derivative_Plot(theta2, omega3, alpha3, dt)   % verify derivatives
Derivative_Plot(theta2, omega4, alpha4, dt)   % verify derivatives
Derivative_Plot(theta2, vP(1,:), aP(1,:), dt) % verify derivatives
 
% Plot the acceleration of point P
figure
plot(theta2*180/pi,aP(1,:),'Color',[153/255 153/255 153/255])
hold on
plot(theta2*180/pi,aP(2,:),'Color',[0/255 153/255 76/255])
legend('a_P_x','a_P_y','Location','Southeast')
title('Acceleration of point P on Fourbar Linkage')
xlabel('Crank angle (\circ)')
ylabel('Acceleration (m/s^2)')
grid on
set(gca,'xtick',0:60:360)
xlim([0 360])

% Save the plot
saveas(gcf, 'Fourbar_Acceleration_Analysis - plot.png')