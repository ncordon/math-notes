
function error = ejercicio2Error(N, a, b, f, truefn)
  approx = ejercicio2Approx(N, a, b, f, truefn);
  h = (b - a)/(N + 1);
  % grid points x including boundaries
  x = linspace(a, b, N + 2);
  % grid points y including boundaries 
  y = linspace(a, b, N + 2);

  % 2d arrays of x,y values
  [X, Y] = meshgrid(x,y);
  X = X';
  Y = Y';
  
  % true solution for test problem
  utrue = arrayfun(truefn, X, Y);
  % assuming true solution is known and stored in utrue:
  error = max(max(abs(approx - utrue)));
endfunction
