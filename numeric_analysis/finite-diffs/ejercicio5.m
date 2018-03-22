
% equation parameter
a = -4;
% Bounds for the space interval
ax = -10;
bx = 2;
% Number of inner nodes
N = 100;
% Initial condition
initCond = @(x) 1 - exp(-x.^2);
f = @(t,x) exp(t) - exp(t - (x + 4*t).^2);

ans = ejercicio5Approx(N, a, ax, bx, initCond, f)
