% Conducts a velocity analysis on the geared fivebar linkage and
% plots the velocity of point P
% by Damian Ziolkowski, February 7, 2021

% Prepare Workspace
clear variables; close all; clc;
 
% Linkage dimensions
a = 0.120;           % crank length (m)
b = 0.250;           % coupler 1 length (m)
c = 0.250;           % coupler 2 length (m)
d = 0.180;           % distance between ground pins (m)
u = 0.120;           % length of link on gear 2 (m)
N1 = 48;             % number of teeth on gear 1
N2 = 24;             % number of teeth on gear 2
rho = N1/N2;         % gear ratio
phi = 0;             % offset angle between gears
gamma3 =  20*pi/180; % angle to point P on coupler 1 (CCW rotation)
gamma4 = -20*pi/180; % angle to point Q on coupler 2 (CW rotation)
p = 0.200;           % distance to point P on coupler 1
q = 0.200;           % distance to point Q on coupler 2

% Ground pins
x0 = [0;0];  % ground pin at A (origin)
xD = [d;0];  % ground pin at D
v0 = [0;0];    % velocity of pin A (zero)

% Angular velocity and acceleration of crank
omega2 = 10;      % angular velocity of crank (rad/s)
omega5 = -rho*omega2;
 
N = 361;   % number of times to perform position calculations
[xB,xC,xE_l,xE_r,xP,xQ] = deal(zeros(2,N));  % allocate space for points
[vB,vC,vE_l,vE_r,vP,vQ] = deal(zeros(2,N));  % allocate space for velocities
[theta2,theta3,theta4]  = deal(zeros(1,N));  % allocate space for angles
                       % "deal" distributes inputs to outputs
[omega3, omega4] = deal(zeros(1,N));  % allocate space for velocities

% Perform calculations for every angle
for i = 1:N
  theta2(i) = (i-1)*(2*pi)/(N-1);  % crank angle
  theta5 = -rho*theta2(i) + phi;   % angle of second gear
  [e2,n2] = UnitVector(theta2(i)); % unit vector for crank
  [e5,n5] = UnitVector(theta5);    % unit vector for second gear
  
  xC(:,i) = FindPos(xD, u, e5);    % coords of pin C 
  
  % From now this is fourbar linkage (ACEB)
  dprime = sqrt(xC(1,i)^2 + xC(2,i)^2);  % distance to pin C (ground pins)
  beta = atan2(xC(2,i),xC(1,i));         % angle to pin C (offset)
  r = dprime - a*cos(theta2(i) - beta);
  s = a*sin(theta2(i) - beta);
  f2 = r^2 + s^2;                        % f squared
  delta = acos((b^2+c^2-f2)/(2*b*c));
  g = b - c*cos(delta);
  h = c*sin(delta);
 
  % Calculate rest of fivebar linkage angles 
  theta3(i) = atan2((h*r - g*s),(g*r + h*s)) + beta;
  theta4(i) = theta3(i) + delta;
  
  % Calculate unit vectors
  [e3,n3] = UnitVector(theta3(i));  % unit vector for first coupler
  [e4,n4] = UnitVector(theta4(i));  % unit vector for second coupler
  [eBP,nBP] = UnitVector(theta3(i) + gamma3); % unit vec from B to P
  [eCQ,nCQ] = UnitVector(theta4(i) + gamma4); % unit vec from C to Q
  
  % Solve for positions of points B, E, P and Q on the linkage
  xB(:,i) = FindPos(     x0,   a,   e2);
   % Point E from left and right arm respectively
  xE_l(:,i) = FindPos(xB(:,i),   b,   e3);
  xE_r(:,i) = FindPos(xC(:,i),   c,   e4);
  xP(:,i) = FindPos(xB(:,i),   p,  eBP);
  xQ(:,i) = FindPos(xC(:,i),   q,  eCQ);    
  
 % Conduct velocity analysis to solve for omega3 and omega4
  A_Mat = [b*n3    -c*n4];
  b_Vec = -a*omega2*n2 + u*omega5*n5;
  omega_Vec = A_Mat\b_Vec;  % solve for angular velocities
 
  omega3(i) = omega_Vec(1);    % decompose omega_Vec into
  omega4(i) = omega_Vec(2);    % individual components
 
% Calculate velocity at important points on linkage
    vB(:,i) = FindVel(     v0,   a,    omega2,   n2);
    vC(:,i) = FindVel(     v0,   u,    omega5,   n5);
  vE_l(:,i) = FindVel(vB(:,i),   b, omega3(i),   n3);
    vP(:,i) = FindVel(vB(:,i),   p, omega3(i),  nBP);
end 

% Code verification - numerical approximation
dt = 2*pi/((N-1)*omega2);             % time increment between calculations
Derivative_Plot(theta2, theta3, omega3, dt)   % verify derivatives
Derivative_Plot(theta2, theta4, omega4, dt)   % verify derivatives
Derivative_Plot(theta2, xP(1,:), vP(1,:), dt) % verify derivatives
 
% Plot the velocity of point P
figure
plot(theta2*180/pi,vP(1,:),'Color',[153/255 153/255 153/255])
hold on
plot(theta2*180/pi,vP(2,:),'Color',[0/255 153/255 76/255])
legend('v_P_x','v_P_y','Location','Southeast')
title('Velocity of point P on Fivebar Linkage')
xlabel('Crank angle (\circ)')
ylabel('Velocity (m/s)')
grid on
set(gca,'xtick',0:60:360)
xlim([0 360])


% Save the plot
saveas(gcf, 'Fivebar_Velocity_Analysis - plot.png')