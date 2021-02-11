% Function Derivative_Plot.m
% plots an approximation to the derivative of a function along with
% the function itself
%
% crankAngle = used for the x axis of the derivative plots
% theta      = the function whose derivative is to be estimated
% omega      = the calculated derivative
% dt         = time step
% sym        = symbol to be printed

function Derivative_Plot(crankAngle, theta, omega, dt)

N = length(theta);       % length of position and velocity vectors
omegaStar = zeros(N,1);  % estimate of derivative

% Estimate derivative
for i = 1:N-1
  omegaStar(i) = (theta(i+1) - theta(i))/dt;
end

% Assume final derivative is the same as the first
omegaStar(N) = omegaStar(1);

% Plot estimated and analytical derivatives
figure
plot(180*crankAngle/pi, omegaStar,'Color', [0/255 153/255 76/255])
hold on
plot(180*crankAngle/pi, omega, '.','Color', [51/255 51/255 51/255])
legend('Estimated','Analytical')
title('Comparison of Calculated and Analytical Derivatives')
xlabel('Crank angle (\circ)')
ylabel('Derivative')
set(gca,'xtick',0:60:360)
xlim([0 360])
grid on

