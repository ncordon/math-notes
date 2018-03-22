
f = @(x, y) 2 * pi * pi * cos(pi * y) * sin(pi * x);
truefn = @(x, y) -cos(pi * y) * sin(pi * x);
% number of points to use
N = 10;
a = 0;
b = 1;

x = linspace(a, b, N + 2);
% grid points y including boundaries 
y = linspace(a, b, N + 2);

% 2d arrays of x,y values
[X, Y] = meshgrid(x, y);
X = X';
Y = Y';

% Compute approximated solution and true solution
approx = ejercicio2Approx(N, a, b, f, truefn);
utrue = arrayfun(truefn, X, Y);

% plot approximation
figure(1, "visible", "off" );
clf;

% plot surface of the original and the approximation values
mesh(X, Y, approx);
title('Surface of approximated solution');
print -dpng ejercicio2-meshapprox.png;
ans = "ejercicio2-meshapprox.png";

% plot true solution
figure(1, "visible", "off" );

% plot surface of the original and the approximation values
mesh(X, Y, utrue);
title('Surface of real solution');
print -dpng ejercicio2-meshsol.png;
ans = "ejercicio2-meshsol.png";
