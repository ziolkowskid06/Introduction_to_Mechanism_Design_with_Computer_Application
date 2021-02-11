% Calculates the unit vector and unit normal for a given angle.
%
% theta = angle of unit vector
% e     = unit vector in the direction of theta
% n     = unit normal to the vector e
 
function [e,n] = UnitVector(theta)
 
e = [ cos(theta); sin(theta)];
n = [-sin(theta); cos(theta)];