
% Approximate solution of
% y''(x) + p(x) y'(x) + q(x)y(x) = f(x)
% y(0) = alpha, y(l) = beta
% using B-splines

N = 30;
ax = 0;
bx = 10;
truesol = @(x) exp(x) .+ x;
alpha = truesol(ax);
beta = truesol(bx);
x = linspace(ax, bx, N + 2);
utrue = @arrayfun(truesol,x);
p = @(x) -2 .+ x .* 0;
q = @(x) 1 .+ x .* 0;
f = @(x) -2 .+ x
approx = ContourApproxEjercicio3(N, ax, bx, p, q, f, alpha, beta);

clf
hold on
plot(x,utrue, x, approx)
legend("true solution", "approximation")
hold off

print -dpng ejercicio3.png;
ans = "ejercicio3.png";
