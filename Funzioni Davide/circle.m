function circle = circle(x,y,r)
% Creates a circle
% Input:
% - x: x coordinate of the center
% - y: y coordinate of the center
% - r: radius of the circle
hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
circle = plot(xunit, yunit);