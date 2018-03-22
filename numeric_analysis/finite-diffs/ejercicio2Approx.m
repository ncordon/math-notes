
% Solve the Poisson problem u_{xx} + u_{yy} = f(x,y) on [a,b] x [a,b]
%
% Based on  http://www.amath.washington.edu/~rjl/fdmbook/chapter3
function approx = ejercicio2Approx(N, a, b, f, truefn)
  h = (b - a)/(N + 1);
  % grid points x including boundaries
  x = linspace(a, b, N + 2);
  % grid points y including boundaries 
  y = linspace(a, b, N + 2);

  % 2d arrays of x,y values
  [X, Y] = meshgrid(x,y);
  X = X';
  Y = Y';

  % indices of interior points in x
  Iint = 2:(N + 1);
  % indices of interior points in y
  Jint = 2:(N + 1);
  % interior points
  Xint = X(Iint, Jint);
  Yint = Y(Iint, Jint);

  % evaluate f at interior points for right hand side
  % rhs is modified below for boundary conditions.
  rhs = arrayfun(f, Xint, Yint);        
  
  % true solution for test problem
  utrue = arrayfun(truefn, X, Y);

  % set boundary conditions around edges of usoln array:
  % use true solution for this test problem
  % This sets full array, but only boundary values
  % are used below.  For a problem where utrue
  % is not known, would have to set each edge of
  % usoln to the desired Dirichlet boundary values.
  usoln = utrue;


  % adjust the rhs to include boundary terms:
  rhs(:, 1) = rhs(:, 1) - usoln(Iint, 1)/h^2;
  rhs(:, N) = rhs(:, N) - usoln(Iint, N + 2)/h^2;
  rhs(1, :) = rhs(1, :) - usoln(1, Jint)/h^2;
  rhs(N, :) = rhs(N, :) - usoln(N + 2, Jint)/h^2;


  % convert the 2d grid function rhs into a column vector for rhs of system:
  % reshape runs through the matrix by columns
  F = reshape(rhs, N * N, 1);

  % form matrix A:
  I = speye(N); % sparse identity matrix
  e = ones(N, 1);
  % [-1 0 1] indicates order of diags, where -1 is the lower-diagonal, 
  % 0 is the main diagonal and 1 is the upper-diagonal
  T = spdiags([e -4*e e], [-1 0 1], N, N);
  S = spdiags([e e], [-1 1], N, N);
  A = (kron(I, T) + kron(S, I)) / h^2;


  % Solve the linear system:
  uvec = A\F;  

  % reshape vector solution uvec as a grid function and 
  % insert this interior solution into usoln for plotting purposes:
  % (recall boundary conditions in usoln are already set) 
  usoln(Iint, Jint) = reshape(uvec, N, N);

  approx = usoln;
endfunction
