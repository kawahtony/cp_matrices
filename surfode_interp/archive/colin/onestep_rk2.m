function vnew = onestep_rk2(f, v, t, k)
%ONESTEP_RK2  Take a single step of the explicit trapezoidal rule
  a = k*f(t, v);
  b = k*f(t + k, v + a);
  vnew = v + 1/2* (a + b);
end
