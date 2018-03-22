function approx = ContourApproxEjercicio4(ax, bx, N, T, M, g)
  % space and time steps
  h = (bx - ax)/(N + 1);
  k = T / M;
  x = linspace(ax, bx, N + 2)';
  t = linspace(0, T, M + 1)';
  % Evaluate g only in inner nodes
  gs = g(x(2:N+1));
  ones_v = ones(1, N);
  ones_v = ones_v';
  
  % Adjust matrixes for the method
  A = spdiags([1/6*ones_v 2/3*ones_v 1/6*ones_v], [-1 0 1], N, N); 
  B = 1/h^2 * spdiags([1*ones_v -2*ones_v 1*ones_v], [-1 0 1], N, N);
  
  a = inv(A) * gs;
  a = a';
  
  for j = 0:M
    % Insert first and last term
    a = [1/3 * a(1), a, 1/3 * a(N)];
    for i = 0:N+1
      bspline = @(q) CubicBspline(q, i, N); 
      approx(i + 1, j + 1) = sum(a .* arrayfun(bspline, 1:N+2));
    endfor
    
    a = a(2:N+1);
    a = (A - k/2 * B) \ (A + k/2 * B) * a';
    a = a';
  endfor
end