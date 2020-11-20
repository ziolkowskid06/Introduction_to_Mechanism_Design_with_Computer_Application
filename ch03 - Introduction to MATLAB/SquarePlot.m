%% Plot a filled square.

clear variables; close all; clc;

% Coordinates of the vertices. First row x-values, second row y-values.
UnitSquare = 0.5*[-1  1 1 -1 -1;
                  -1 -1 1  1 -1];

% Enlarge square dimensions twice and center the square at (1, 1).
NewSquare = 2*UnitSquare + [1; 1];

% Fill area inside the square.
fill(NewSquare(1,:), NewSquare(2,:), 'g', 'LineWidth', 2)

% Plot parameters.
title('Filled square')
xlabel('x (m)'); xlim([-2 2])
ylabel('y (m)'); ylim([-2 2])
grid on
axis equal

% Save the plot.
saveas(gcf, 'SquarePlot.png')