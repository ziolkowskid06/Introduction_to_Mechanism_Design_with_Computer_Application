% Function FindVel.m
% calculates the translational velocity at a point on the linkage
% using the relative velocity formula
%
% v0    = velocity of first point
% L     = length of vector to second point on the link
% omega = angular velocity of link
% n     = unit normal to vector btw first and second points
% v     = velocity of second point
 
function v = FindVel(v0, L, omega, n)
 
v = v0 + omega * L * n;