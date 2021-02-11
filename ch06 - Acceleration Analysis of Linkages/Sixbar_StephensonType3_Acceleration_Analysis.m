% Conducts an acceleration analysis of the Stephenson Type III sixbar linkage
% and plots the acceleration of point E.
% by Damian Ziolkowski, February 10, 2021

% Prepare Workspace
clear variables; close all; clc;
 
% Linkage dimensions
a = 0.080;           % AB crank length (m)
b = 0.105;           % BC coupler length (m)
c = 0.120;           % CD rocker length (m)
d = 0.135;           % AD length between ground pins (m)
p = 0.225;           % BE length on coupler (m)
q = 0.195;           % AF length on ground (m)
u = 0.180;           % EG length of link 5 (m)
v = 0.180;           % FG length of link 6 (m)
gamma1 = -20*pi/180; % internal angle of ground (CW rotation)
gamma3 = 20*pi/180;  % internal angle of coupler triangle (CCW rotation)

% Ground pins
x0 = [ 0; 0];    % ground pin at A (origin)
xD = [ d; 0];    % ground pin at D
[eAF,nAF] = UnitVector(gamma1);  
xF = FindPos(x0, q, eAF);   % ground pin at F
v0 = [0; 0];     % velocity of pin A (zero)
a0 = [0; 0];     % acceleration of pin A (zero)

% Angular velocity and acceleration of crank
omega2 = 10;      % angular velocity of crank (rad/s)
alpha2 = 0;       % angular acceleration of crank (rad/s^2)

% Allocate space for variables
N = 361;   % number of times to perform position calculations
Z21 = zeros(2,1);   % column vector of two zeros
[xB,xC,xE,xG]                        = deal(zeros(2,N)); % points
[vB,vC,vE,vG]                        = deal(zeros(2,N)); % points
[aB,aC,aE,aG]                        = deal(zeros(2,N)); % acceleration
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
  [eBE,nBE] = UnitVector(theta3(i) + gamma3);
  
  % Solve for positions of points B, C, E
  xB(:,i) = FindPos(x0, a, e2);
  xC(:,i) = FindPos(xB(:,i), b, e3);
  xE(:,i) = FindPos(xB(:,i), p, eBE);
  
  % Solve right fourbar linkage
  xFC = xF(1) - xC(1,i);   
  yFC = xF(2) - xC(2,i);  
  xEC = xE(1,i) - xC(1,i);   
  yEC = xE(2,i) - xC(2,i);
  beta =  atan2(yFC, xFC);
  alpha = atan2(yEC, xEC);
  theta2Prime = alpha - beta;          % virtual crank angle on upper fourbar
  aPrime = sqrt(xEC^2 + yEC^2);        % virtual crank length on right fourbar
  dPrime = sqrt(xFC^2 + yFC^2);        % virtual ground length on upper fourbar
 
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
  A_Mat = [b*n3    -c*n4    Z21    Z21; p*nBE  Z21   u*n5   -v*n6];
  b_Vec = [-a*omega2*n2; -a*omega2*n2];
  omega_Vec = A_Mat\b_Vec;  % solve for angular velocities
 
  omega3(i) = omega_Vec(1);    % decompose omega_Vec into
  omega4(i) = omega_Vec(2);    % individual components
  omega5(i) = omega_Vec(3);   
  omega6(i) = omega_Vec(4);
 
% Calculate velocity at important points on linkage
  vB(:,i) = FindVel(     v0,   a,    omega2,   n2);
  vC(:,i) = FindVel(     v0,   c, omega4(i),   n4);
  vE(:,i) = FindVel(vB(:,i),   p, omega3(i),  nBE);
  vG(:,i) = FindVel(     v0,   v, omega6(i),   n6);
    
% Conduct acceleration analysis to solve for alpha3, alpha4, alpha5, alpha6
  at = a*alpha2;         % tangential acceleration  
  ac = a*omega2^2;       % centripetal acceleration 
  bc = b*omega3(i)^2;    % centripetal acceleration 
  cc = c*omega4(i)^2;    % centripetal acceleration 
  pc = p*omega3(i)^2;    % centrupetal acceleration
  uc = u*omega5(i)^2;    % centripetal acceleration 
  vc = v*omega6(i)^2;    % centripetal acceleration
 
  C_Mat = A_Mat;
  d_Vec = [-at*n2 + ac*e2 + bc*e3 - cc*e4; -at*n2 + ac*e2 + pc*eBE + uc*e5 - vc*e6];
  alpha_Vec = C_Mat\d_Vec;  % solve for angular accelerations
 
  alpha3(i) = alpha_Vec(1);
  alpha4(i) = alpha_Vec(2);
  alpha5(i) = alpha_Vec(3);
  alpha6(i) = alpha_Vec(4);
 
% Find acceleration of pins
  aB(:,i) = FindAcc(     a0, a,    omega2,    alpha2,  e2,  n2);
  aC(:,i) = FindAcc(     a0, c, omega4(i), alpha4(i),  e4,  n4);
  aE(:,i) = FindAcc(aB(:,i), p, omega3(i), alpha3(i), eBE, nBE);
  aG(:,i) = FindAcc(aE(:,i), u, omega5(i), alpha5(i),  e5,  n5);
end

% Code verification - numerical approximation
dt = 2*pi/((N-1)*omega2);             % time increment between calculations
Derivative_Plot(theta2, theta3, omega3, dt)   % verify derivatives
Derivative_Plot(theta2, theta4, omega4, dt)   % verify derivatives
Derivative_Plot(theta2, theta5, omega5, dt)   % verify derivatives
Derivative_Plot(theta2, theta6, omega6, dt)   % verify derivatives
Derivative_Plot(theta2,vE(1,:),aE(1,:), dt)   % verify derivatives
 
% Plot the acceleration of point E
figure
plot(theta2*180/pi,aE(1,:),'Color',[153/255 153/255 153/255])
hold on
plot(theta2*180/pi,aE(2,:),'Color',[0/255 153/255 76/255])
legend('a_E_x','a_E_y','Location','Southeast')
title('Acceleration of point E on Stephenson Type III Linkage')
xlabel('Crank angle (\circ)')
ylabel('Acceleration (m/s^2)')
grid on
set(gca,'xtick',0:60:360)
xlim([0 360])

% Save the plot
saveas(gcf, 'Sixbar_StephensonType3_Acceleration_Analysis - plot.png')
