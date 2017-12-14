function [d,x] = radSlope(x,y)
% Calculates slope between 3 points, starting from 2 in the array. 1 is
% zero

n = length(x);
d = zeros (1,n);

dx = x(2)-x(1);
dy = y(2)-y(1);

% x = 1:1/2:n/2+1;
x = linspace(1,n+1,n*8);

x0 = x(1);
y0 = y(1);

xn = x(n);
yn = y(n);

d(1) = atan((y0-y(2))/(x(2)-x0));

for i=2:1:n-1
    d(i) = atan((y(i-1)-y(i+1))/(x(i+1)-x(i-1)));
end

d(n) = atan((y(n-1)-yn)/(xn-x(n-1)));

% d = pi-d;