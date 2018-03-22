
% solve the periodic boundary value problem 
%       u'' - u = f(t) on [0,1] 
% with Dirichlet conditions u(0) = u(1), u'(0) = u'(1).  
% Modified from http://www.amath.washington.edu/~rjl/fdmbook/chapter3

f = @(t) cos(pi * t);
% True solution to Poisson problem
truefn = @(t) (exp(1-t) - exp(t) -(-1+e)*f(t)) / ((e-1)*(1 + pi*pi));

% number of interior points 
N = 20; 
h = 1/(N + 1);
a = -2 - h*h;

% grid points x including boundaries
x = linspace(0, 1, N+1);
X = x';
rhs = f(X);
              
% form matrix A:
A = sparse(N+1, N+1);
% Assign diagonals
A(1:N+2:end) = a;
A(2:N+2:end) = 1;
A(N+2:N+2:end) = 1;
% Assign firs/last 1 of last/first row resp.
A(1 + (N+1)*N) = A(N+1) = 1;
A *= 1/(h*h);

% Approximated solution
uapprox = A\rhs;

% True solution for test problem
utrue = truefn(X);

err = max(max(abs(uapprox - utrue)));
fprintf('Error relative to true solution of PDE = %10.3e \n', err);


% plot results:
figure( 1, "visible", "off" );
clf;
hold on;

% plot true and aproximated solution:
plot(X, utrue, X, uapprox);
title('Plot of approximated and true solution');
legend("True solution", "Approximated solution");
print -dpng ejercicio1.png;
ans = "ejercicio1.png";
