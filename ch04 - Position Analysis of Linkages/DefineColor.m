% Makes a palette of custom colors for plotting in MATLAB
% the first color in the palette gives the full base color, and
% the rest fade gradually to white
%
% colorBase = 1x3 input array giving RGB values (between 0 and 255)
% C         = 10x3 output array giving RGB values (between 0 and 1)
 
function C = DefineColor(colorBase)
 
C = zeros(11,3);
for i = 1:11
  for j = 1:3
    C(i,j) = (i-1)*(255 - colorBase(j))/10 + colorBase(j);
  end
end
 
C = C/255; % MATLAB specifies its RGB values in the range between 0 and 1 