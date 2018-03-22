
function approx = ContourApproxEjercicio3(N, ax, bx, p, q, f, alpha, beta)
  h = (bx - ax)/(N + 1);
  x = linspace(ax, bx, N + 2)';
  
  lower_diag =  zeros(1, N + 2);
  middle_diag = zeros(1, N + 2);
  upper_diag =  zeros(1, N + 2);
  
  qs = q(x)';
  ps = p(x)';
  fs = f(x)';
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Compute middle diagonal
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  pmid_coeffs = zeros(1, N + 2);
  pmid_coeffs(1) = 3/h;
  pmid_coeffs(2) = 1/(4 * h);
  pmid_coeffs(N + 2) = -3/h;
  pmid_coeffs(N + 1) = -1/(4 * h);
  
  middle_diag += pmid_coeffs .* ps;
  
  hhmid_coeffs = ones(1, N+2) * -2/h^2;
  hhmid_coeffs(1) = hhmid_coeffs(N + 2) = -9/h^2;
  hhmid_coeffs(2) = hhmid_coeffs(N + 1) = -5/(2 * h^2);
  
  middle_diag += hhmid_coeffs;
  
  qmid_coeffs = ones(1, N+2) * 2/3;
  qmid_coeffs(1) = qmid_coeffs(N + 2) = 0;
  qmid_coeffs(2) = qmid_coeffs(N + 1) = 7/12;
  
  middle_diag += qmid_coeffs .* qs;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Compute lower diagonal
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  plower_coeffs = ones(1, N + 2) * -1/(2 * h);
  plower_coeffs(2) = -3/(4 * h);
  plower_coeffs(N + 2) = 0;
  plower_coeffs = plower_coeffs .* ps;
  plower_coeffs = [plower_coeffs(2:N+2), 0];
  
  lower_diag += plower_coeffs;
  
  hhlower_coeffs = ones(1, N + 2) * 1/h^2;
  hhlower_coeffs(2) = 3/(2 * h^2);
  hhlower_coeffs(N + 2) = 3/h^2;
  hhlower_coeffs = [hhlower_coeffs(2:N+2), 0];
  
  lower_diag += hhlower_coeffs;
  
  qlower_coeffs = ones(1, N + 2) * 1/6;
  qlower_coeffs(2) = 1/4;
  qlower_coeffs(N + 2) = 0;
  qlower_coeffs = qlower_coeffs .* qs;
  qlower_coeffs = [qlower_coeffs(2:N+2), 0];
  
  lower_diag += qlower_coeffs;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Compute upper diagonal
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  pupper_coeffs = ones(1, N + 2) * 1/(2 * h);
  pupper_coeffs(1) = 0;
  pupper_coeffs(N + 1) = 3/(4 * h);
  pupper_coeffs = pupper_coeffs .* ps;
  pupper_coeffs = [0, pupper_coeffs(1:N+1)];
  
  upper_diag += pupper_coeffs;
  
  hhupper_coeffs = ones(1, N + 2) * 1/h^2;
  hhupper_coeffs(1) = 3/h^2;
  hhupper_coeffs(N + 1) = 3/(2 * h^2);
  hhupper_coeffs = [0, hhupper_coeffs(1:N+1)];
  
  upper_diag += hhupper_coeffs;
  
  qupper_coeffs = ones(1, N + 2) * 1/6;
  qupper_coeffs(1) = 0;
  qupper_coeffs(N + 1) = 1/4;
  qupper_coeffs = qupper_coeffs .* qs;
  qupper_coeffs = [0, qupper_coeffs(1:N+1)];
  
  upper_diag += qupper_coeffs;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Compute coefficients of linear combination of B-splines
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  A = spdiags([lower_diag' middle_diag' upper_diag'], [-1 0 1], N + 2, N + 2);
  rhs = fs;
  rhs(1) += -q(1)*alpha + 3*p(1)*alpha/h - 6/h^2 * alpha;
  rhs(N + 2) += -q(N + 2)*beta - 3*p(N + 2)*beta/h - 6/h^2 * beta;
  
  coeffs = A\rhs';
  coeffs = [alpha, coeffs', beta];
  approx = 0 .* x;
  
  for i = 1:N+4
    bspline = @(t) CubicBspline(i - 1, t/h, N);
    approx += coeffs(i) .* arrayfun(bspline, x);
  endfor
end