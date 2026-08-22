function vnew = onestep_rk4(f, v, t, k)
%ONESTEP_RK4  Take a single step of the classical 4-stage Runge-Kutta method
%   VNEW = ONESTEP_RK4(F, V, T, K)
%         Inputs:
%           F: a function handle that takes (t, u) as inputs
%           V: current state
%           T: current time
%           K: timestep
%         Returns:
%           VNEW: solution at next step.
  a = k*f(t, v);
  b = k*f(t + k/2, v + a/2);
  c = k*f(t + k/2, v + b/2);
  d = k*f(t + k, v + c);
  vnew = v + 1/6* (a + 2*b + 2*c + d);
end
