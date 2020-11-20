%% Makes a plot of eight squares arranged in a circle.

clear variables; close all; clc

UnitSquare = 0.5*[-1  1 1 -1 -1;
                  -1 -1 1  1 -1];

% Prepare eight squares 
for i = 1:8
  theta = (i-1)*2*pi/8;           % angle at which to place square
  r = 2*[cos(theta); sin(theta)]; % center of square
  NewSquare = UnitSquare + r;     % create new square
  col = (i-1)*[1 1 1]/8;          % fill and edge color of square
  fill(NewSquare(1,:),NewSquare(2,:),col,'EdgeColor',col,'LineWidth',2)
  hold on           % Keep all squares on a plot
end

% Plot parameters
title('Eight squares arranged in a circle')
xlabel('x (m)'); xlim([-3 3]) 
ylabel('y (m)'); ylim([-3 3])
grid on
axis equal

% Save plot
saveas(gcf, 'SquareCircle.png')