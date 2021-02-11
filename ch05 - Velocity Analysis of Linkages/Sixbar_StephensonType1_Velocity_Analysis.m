% Conducts a velocity analysis on the Stephenson Type I sixbar linkage
% and plots the velocity of point G.
% by Damian Ziolkowski, February 7, 2021
 
% Prepare Workspace
clear variables; close all; clc;
 
% Linkage dimensions
a = 0.070;           % crank length (m)
b = 0.100;           % coupler length (m)
c = 0.090;           % rocker length (m)
d = 0.110;           % length between ground pins (m)
p = 0.150;           % length to third pin on crank triangle (m)
q = 0.150;           % length to third pin on rocker triangle (m)
u = 0.120;           % length of link 5 (m)
v = 0.160;           % length of link 6 (m)
gamma2 = 20*pi/180;  % internal angle of crank triangle (CCW rotation)
gamma4 = -20*pi/180; % internal angle of rocker triangle (CW rotation)

% Ground pins
x0 = [ 0; 0];    % ground pin at A (origin)
xD = [ d; 0];    % ground pin at D
v0 = [0;0];      % velocity of pin A (zero)

% Angular velocity and acceleration of crank
omega2 = 10;      % angular velocity of crank (rad/s)
 
% Allocate space for variables
N = 361;   % number of times to perform position calculations
Z21 = zeros(2,1);   % column vector of two zeros
[xB,xC,xE,xF,xG]                     = deal(zeros(2,N)); % points
[vB,vC,vE,vF,vG]                     = deal(zeros(2,N)); % velocities
[theta2,theta3,theta4,theta5,theta6] = deal(zeros(1,N)); % angles
                                    % "deal" distributes inputs to outputs
[omega3, omega4, omega5, omega6] = deal(zeros(1,N));  % angular velocities

% Perform calculations for every angle
for i = 1:N
 
  % Solve lower fourbar linkage
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
  [eAE,nAE] = UnitVector(theta2(i) + gamma2);
  [eDF,nDF] = UnitVector(theta4(i) + gamma4);
 
  % Solve for positions of points B, C, E, F
  xB(:,i) = FindPos(x0, a,  e2);
  xC(:,i) = FindPos(xD, c,  e4);
  xE(:,i) = FindPos(x0, p, eAE);
  xF(:,i) = FindPos(xD, q, eDF);
  
  % Solve upper fourbar linkage
  xFB = xF(1,i) - xB(1,i);   yFB = xF(2,i) - xB(2,i);  
  xEB = xE(1,i) - xB(1,i);   yEB = xE(2,i) - xB(2,i);
  beta =  atan2(yFB, xFB);
  alpha = atan2(yEB, xEB);
  aPrime = sqrt(xEB^2 + yEB^2);        % virtual crank length on upper fourbar
  dPrime = sqrt(xFB^2 + yFB^2);        % virtual ground length on upper fourbar
  theta2Prime = alpha - beta;          % virtual crank angle on upper fourbar
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
  
  % Conduct velocity analysis to solve for omega3, omega4, omega5, omega6
  A_Mat = [b*n3    -c*n4    Z21    Z21; Z21   -q*nDF   u*n5   -v*n6];
  b_Vec = [-a*omega2*n2; -p*omega2*nAE];
  omega_Vec = A_Mat\b_Vec;  % solve for angular velocities
 
  omega3(i) = omega_Vec(1);    % decompose omega_Vec into
  omega4(i) = omega_Vec(2);    % individual components
  omega5(i) = omega_Vec(3);   
  omega6(i) = omega_Vec(4);
 
% Calculate velocity at important points on linkage
  vB(:,i) = FindVel(     v0,   a,    omega2,   n2);
  vC(:,i) = FindVel(     v0,   c, omega4(i),   n4);
  vF(:,i) = FindVel(     v0,   q, omega4(i),  nDF);
  vE(:,i) = FindVel(     v0,   p,    omega2,  nAE);
  vG(:,i) = FindVel(vE(:,i),   u, omega5(i),   n5);
end

% Code verification - numerical approximation
dt = 2*pi/((N-1)*omega2);             % time increment between calculations
Derivative_Plot(theta2, theta3, omega3, dt)   % verify derivatives
Derivative_Plot(theta2, theta4, omega4, dt)   % verify derivatives
Derivative_Plot(theta2, theta5, omega5, dt)   % verify derivatives
Derivative_Plot(theta2, theta6, omega6, dt)   % verify derivatives
Derivative_Plot(theta2,xG(1,:),vG(1,:), dt)   % verify derivatives
 
% Plot the velocity of point G
figure
plot(theta2*180/pi,vG(1,:),'Color',[153/255 153/255 153/255])
hold on
plot(theta2*180/pi,vG(2,:),'Color',[0/255 153/255 76/255])
legend('v_G_x','v_G_y','Location','Southeast')
title('Velocity of point G on Stephenson Type I Linkage')
xlabel('Crank angle (\circ)')
ylabel('Velocity (m/s)')
grid on
set(gca,'xtick',0:60:360)
xlim([0 360])

% Save the plot
saveas(gcf, 'Sixbar_StephensonType1_Velocity_Analysis - plot.png')
