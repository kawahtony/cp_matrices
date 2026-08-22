function vnew = onestep_heun(f, v, t, k)
%ONESTEP_HEUN  Take a single step of Heun's method
  a = k*f(t, v);
  b = k*f(t + k/3, v + a/3);
  c = k*f(t + 2*k/3, v + 2*b/3);
  vnew = v + 1/4*(a + 3*c);
end
