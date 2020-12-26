% Calculates the position of a point on a link using the ...
% relative position formula
%
% x0 = position of first point on the link
% L  = length of vector between first and second points
% e  = unit vector between first and second points
% x  = position of second point on the link
 
function x = FindPos(x0, L, e)
 
x = x0 + L * e;