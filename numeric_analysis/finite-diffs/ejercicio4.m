
function [x, uint] = ejercicio4(N, a, ax, bx, initCond, f)
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
  % Frequency of plot
  nplot = 20;
  
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
    u(I) = u(I) - a * k / (2 * h) * (3*u(I) - 4*u(I - 1) + U(I - 2)) + ...
           a^2 * k^2 / (2 * h^2) * (u(I) - 2*u(I - 1) + U(I - 2));
    
    % plot results each nplot steps
    if mod(n, nplot) == 0 || n == nsteps
      uint = u(1:m+2);
      plot(x, uint)
      axis([0 1 -.2 1.2])
      title(sprintf('t = %9.5e  tras %4i pasos con %5i nodos',...
                    tnp, n, m+1))
      if n < nsteps, pause(0.2); end;
    end

    % Periodic conditions
    u(1) = u(N + 2);   
    u(N + 3) = u(2);

    % Set next temporal step
    tn = tnp;
  endfor
  
endfunction
