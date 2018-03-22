
% Approximate solution of
%
% u_{xx} = u_{t}, 0 < x < l; t > 0
% u(0,t) = 0, t >= 0;  u(l,t) = 0, t >= 0
% u(x,0) = g(x), 0 <= x <= l
%
% using B-splines

% Number of space and time evolutions
N = 10;
M = 100;
% Time and space intervals
ax = 0;
bx = 1;
T = 2;
g = @(x) sin(pi * x);
truesol = @(x,t) exp(-pi^2 * t) .* sin(pi * x);
x = linspace(ax, bx, N + 2)';
t = linspace(0, T, M + 1)';
[xs, ts] = meshgrid(x, t);
utrue = truesol(xs, ts);
utrue = utrue';

approx = ContourApproxEjercicio4(ax, bx, N, T, M, g);
%mesh(xs', ts', utrue);
mesh(xs', ts', approx);
title("Approximation")
print -dpng ejercicio4.png;
ans = "ejercicio4.png";
