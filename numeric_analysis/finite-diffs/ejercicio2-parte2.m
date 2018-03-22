
% plot approximation errors

f = @(x, y) 2 * pi * pi * cos(pi * y) * sin(pi * x);
truefn = @(x, y) -cos(pi * y) * sin(pi * x);
% number of points to use
a = 0;
b = 1;

Ns = 1:1:100;
hs = 1./(Ns + 1);
errors = arrayfun(@ejercicio2Error, Ns, a, b, f, truefn);

figure(1, "visible", "off" );
clf;
hold on;
plot(hs, errors, 'k');
xlabel('h value');
ylabel('maximum error');
title('Approximation errors');
print -dpng ejercicio2.png;
ans = "ejercicio2.png";
