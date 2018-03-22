% 0 <= t <= N+1
function value = CubicBspline(i, t, N)
  if i == 0
    value = TPower(1, t, 3);
  elseif i == 1
    value = 1/4 * TPower(2, t, 3) - 2 * TPower(1, t, 3);
  elseif i == 2
    value = 1/6 * TPower(3, t, 3) - 3/4 * TPower(2, t, 3) + 3/2 * TPower(1, t, 3);
  elseif i == 3
    value = 1/6 * TPower(4, t, 3) - 4/6 * TPower(3, t, 3) + TPower(2, t, 3) - 4/6 * TPower(1, t, 3);
  elseif i > 3 && i <= N
    value = CubicBspline(3, t - (i-3), N);
  elseif i > N && i <= N + 3
    value = CubicBspline((N + 3) - i, N + 1 - t, N);
  endif
endfunction