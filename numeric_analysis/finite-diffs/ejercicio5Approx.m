
function ans = ejercicio5Approx(N, a, ax, bx, initCond, f)
  % maximum time
  tfinal = 1;

  % Discretizations
  h = (bx - ax)/(N + 1);
  k = 0.05*h;
  % Courant parameter
  nu = a*k/h;
  % Recall that x(1) = 0, x(N + 2) = 1
  x = linspace(ax, bx, N + 2)';
  % Number of time steps
  nsteps = round(tfinal / k);

  % With periodic conditions we have N + 1 unknowns and u(1) = u(N + 2)
  I = 2:(N + 2); 

  % Definition of initial conditions
  tn = 0;
  u0 = initCond(x);
  u = u0;
  xs = x;
  xs(N + 3) = bx + k;
  
  % Periodic boundary conditions
  u(1) = u(N + 2);   
  u(N + 3) = u(2);

  % Plot initial condition
  clf
  plot(x, u0)
  axis([0 1 -.2 1.2])
  title('Initial condition at 0 time')
  
  % Temporal loop
  for n = 1:nsteps  
    tnp = tn + k;   % = t_{n+1}
    % Compute values for f_i^j
    fs = f(tn, xs);
    nextfs = f(tnp, xs);
    
    % Lax-Wendroff
    u(I) = u(I) - 0.5*nu*(u(I + 1) - u(I - 1)) + ...
           0.5*nu^2 * (u(I - 1) - 2*u(I) + u(I + 1)) - ...
           a * k^2 / (4 * h) * (fs(I + 1) - fs(I - 1)) + ...
           0.5* k * (nextfs(I) - fs(I)) + k * fs(I);
    
    % Periodic conditions
    u(1) = u(N + 2);   
    u(N + 3) = u(2);

    % Set next temporal step
    tn = tnp;
  endfor
  
  uint = u(1:N+2);
  plot(x, uint)
  title(sprintf('t = %9.5e  after %4i steps, %5i nodes', tnp, nsteps, N + 1))
  print -dpng ejercicio5.png;
  ans = "ejercicio5.png";
endfunction
