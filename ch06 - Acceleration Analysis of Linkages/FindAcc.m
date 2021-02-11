% Function FindAcc.m
% Calculates the translational acceleration of a point on the linkage
% using the relative acceleration formula
%
% a0    = acceleration of first point
% L     = length of vector to second point on the link
% omega = angular velocity of link
% alpha = angular acceleration of link
% e     = unit vector btw first and second points
% n     = unit normal to vector btw first and second points
% a     = acceleration of second point
 
function a = FindAcc(a0, L, omega, alpha, e, n)
 
a = a0 + L*alpha*n - L*omega^2*e;