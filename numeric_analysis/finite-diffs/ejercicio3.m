
% Solve u_t = u_{xx} on [ax,bx] with Dirichlet boundary conditions
%
% Based on http://www.amath.washington.edu/~rjl/fdmbook/heat_CN.m



% Number of interior points
N = 20;
% [0,1] space interval
ax = 0;           
bx = 1;
% [0,2] time interval
tfinal = 2;

h = (bx - ax)/(N + 1);
x = linspace(ax, bx, N + 2)';

% time step, 4* not necessary, since there's no stability condition
k = 4*h;
% number of time steps
nsteps = round(tfinal / k);


% true solution for comparison and definition of the boundary conditions
utrue = @(t, x) exp(-pi * pi * t) * sin(pi * x);

% initial conditions:
u0 = utrue(0, x);


% Each time step we solve MOL system U' = AU + g using the Trapezoidal method

% set up matrixes:
r = k/(h^2);
e = ones(N, 1);
A = spdiags([e -2*e e], [-1 0 1], N, N);
A1 = eye(N) - r * A;
A2 = eye(N);


% initial data on fine grid for plotting:
xfine = linspace(ax, bx, 1001);
ufine = utrue(0, xfine);

% initialize u and plot:
tn = 0;
u = u0;

clf;
hold on;

plot(x, u, 'b.-', xfine, ufine, 'r');
legend('computed', 'true');

% main time-stepping loop:
for n = 1:nsteps
  tnp = tn + k;   % = t_{n+1}

  % boundary values u(0,t) = g_0(t) and u(1,t) = g_1(t) at times tn and tnp:
  g0n = u(1);
  g1n = u(N + 2);
  g0np = utrue(tnp, ax);
  g1np = utrue(tnp, bx);

  % compute right hand side for linear system:
  uint = u(2:(N + 1));   % interior points (unknowns)
  rhs = A2*uint;
  % fix-up right hand side using BC's (i.e. add vector g to A2*uint)
  rhs(1) = rhs(1) + r*(g0n + g0np);
  rhs(N) = rhs(N) + r*(g1n + g1np);

  % solve linear system:
  uint = A1\rhs;

  % augment with boundary values:
  u = [g0np; uint; g1np];

  %Plot results
  ufine = utrue(tnp, xfine);
  plot(x, u, 'b.-', xfine, ufine, 'r');

  % for next time step
  tn = tnp;
end % del for


print -dpng ejercicio3.png;
ans = "ejercicio3.png";
