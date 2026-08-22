function vnew = onestep_euler(f, v, t, k)
%ONESTEP_EULER  Take a single step of forward Euler
  a = k*f(t, v);
  vnew = v + a;
end
