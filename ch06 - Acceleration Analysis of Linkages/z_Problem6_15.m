% Problem 6.15 (Watt Type II sixbar linkage)
% by Damian Ziolkowski, February 11, 2021

% Prepare Workspace
clear variables; close all; clc;
 
% Linkage dimensions
a = 0.070;           % AB crank length (m)
b = 0.100;           % BC coupler length (m)
c = 0.090;           % CD rocker length (m)
d = 0.110;           % AD length between ground pins (m)
p = 0.150;           % DE length on coupler (m)
q = 0.150;           % AF length on rocker (m)
u = 0.120;           % EG length of link 5 (m)
v = 0.160;           % FG length of link 6 (m)
z = 0.100;           % EP length of link 7 (m)
gamma4 = -30*pi/180; % internal angle of rocker (CW rotation)
gamma3 =  50*pi/180; % internal angle of triangle EGP (CCW rotation)
gamma1 = -20*pi/180; % internal angle of ground (CW rotation)

% Ground pins
x0 = [ 0; 0];    % ground pin at A (origin)
xD = [ d; 0];    % ground pin at D
[eAF,nAF] = UnitVector(gamma1);  
xF = FindPos(x0, q, eAF);   % ground pin at F
v0 = [0;0];      % velocity of pin A (zero)
a0 = [0;0];      % acceleration of pin A (zero)

% Angular velocity and acceleration of crank
omega2 = 10;      % angular velocity of crank (rad/s)
alpha2 = 0;       % angular acceleration of crank (rad/s^2)

% Allocate space for variables
N = 361;   % number of times to perform position calculations
Z21 = zeros(2,1);   % column vector of two zeros
[xB,xC,xE,xG,xP]                     = deal(zeros(2,N)); % points
[vB,vC,vE,vG,vP]                     = deal(zeros(2,N)); % velocities
[aB,aC,aE,aG,aG]                     = deal(zeros(2,N)); % accelerations
[theta2,theta3,theta4,theta5,theta6] = deal(zeros(1,N)); % angles
                                    % "deal" distributes inputs to outputs
[omega3, omega4, omega5, omega6] = deal(zeros(1,N));  % angular velocities
[alpha3, alpha4, alpha5, alpha6] = deal(zeros(1,N));  % angular accelerations

% Perform calculations for every angle
for i = 1:N
 
  % Solve left fourbar linkage
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
  [eDE,nDE] = UnitVector(theta4(i) + gamma4);
 
  % Solve for positions of points B, C, E, F
  xB(:,i) = FindPos(x0, a,  e2);
  xC(:,i) = FindPos(xD, c,  e4);
  xE(:,i) = FindPos(xD, p, eDE);
  
  % Solve right fourbar linkage
  xFD = xF(1) - xD(1);   
  yFD = xF(2) - xD(2);  
  xED = xE(1,i) - xD(1);   
  yED = xE(2,i) - xD(2);
  beta =  atan2(yFD, xFD);
  alpha = atan2(yED, xED);
  theta2Prime = alpha - beta;          % virtual crank angle on upper fourbar
  aPrime = sqrt(xED^2 + yED^2);        % virtual crank length on right fourbar
  dPrime = sqrt(xFD^2 + yFD^2);        % virtual ground length on upper fourbar
 
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
  [eEP,nEP] = UnitVector(theta5(i) + gamma3);
                 
  % Calculate position of point G
  xG(:,i) = FindPos(xE(:,i), u, e5); 
  xP(:,i) = FindPos(xE(:,i), z, eEP);
  
  % Conduct velocity analysis to solve for omega3, omega4, omega5, omega6
  A_Mat = [b*n3    -c*n4    Z21    Z21; Z21   p*nDE   u*n5   -v*n6];
  b_Vec = [-a*omega2*n2; Z21];
  omega_Vec = A_Mat\b_Vec;  % solve for angular velocities
 
  omega3(i) = omega_Vec(1);    % decompose omega_Vec into
  omega4(i) = omega_Vec(2);    % individual components
  omega5(i) = omega_Vec(3);   
  omega6(i) = omega_Vec(4);
 
% Calculate velocity at important points on linkage
  vB(:,i) = FindVel(     v0,   a,    omega2,   n2);
  vC(:,i) = FindVel(     v0,   c, omega4(i),   n4);
  vE(:,i) = FindVel(     v0,   p, omega4(i),  nDE);
  vG(:,i) = FindVel(     v0,   v, omega6(i),   n6);
  vP(:,i) = FindVel(vE(:,i),   z, omega5(i),  nEP);
  
% Conduct acceleration analysis to solve for alpha3, alpha4, alpha5, alpha6
  at = a*alpha2;         % tangential acceleration  
  ac = a*omega2^2;       % centripetal acceleration 
  bc = b*omega3(i)^2;    % centripetal acceleration 
  cc = c*omega4(i)^2;    % centripetal acceleration 
  pc = p*omega4(i)^2;    % centripetal acceleration
  uc = u*omega5(i)^2;    % centripetal acceleration 
  vc = v*omega6(i)^2;    % centripetal acceleration
 
  C_Mat = A_Mat;
  d_Vec = [-at*n2 + ac*e2 + bc*e3 - cc*e4; pc*eDE + uc*e5 - vc*e6];
  alpha_Vec = C_Mat\d_Vec;  % solve for angular accelerations
 
  alpha3(i) = alpha_Vec(1);
  alpha4(i) = alpha_Vec(2);
  alpha5(i) = alpha_Vec(3);
  alpha6(i) = alpha_Vec(4);
 
% Find acceleration of pins
  aB(:,i) = FindAcc(     a0, a,    omega2,    alpha2,  e2,  n2);
  aC(:,i) = FindAcc(     a0, c, omega4(i), alpha4(i),  e4,  n4);
  aE(:,i) = FindAcc(     a0, p, omega4(i), alpha4(i), eDE, nDE);
  aG(:,i) = FindAcc(aE(:,i), u, omega5(i), alpha5(i),  e5,  n5);
  aP(:,i) = FindAcc(aE(:,i), z, omega5(i), alpha5(i), eEP, nEP);
end

% Code verification - numerical approximation
dt = 2*pi/((N-1)*omega2);             % time increment between calculations
Derivative_Plot(theta2, omega3, alpha3, dt)   % verify derivatives
Derivative_Plot(theta2, omega4, alpha4, dt)   % verify derivatives
Derivative_Plot(theta2, omega5, alpha5, dt)   % verify derivatives
Derivative_Plot(theta2, omega6, alpha6, dt)   % verify derivatives
Derivative_Plot(theta2,vP(1,:),aP(1,:), dt)   % verify derivatives
 
% Plot the acceleration of point P
figure
plot(theta2*180/pi,aP(1,:),'Color',[153/255 153/255 153/255])
hold on
plot(theta2*180/pi,aP(2,:),'Color',[0/255 153/255 76/255])
legend('a_P_x','a_P_y','Location','Northeast')
title('Acceleration of point P on Watt Type II Linkage')
xlabel('Crank angle (\circ)')
ylabel('Acceleration (m/s^2)')
grid on
set(gca,'xtick',0:60:360)
xlim([0 360])