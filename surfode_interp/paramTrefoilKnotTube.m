function [x,y,z,u,v] = paramTrefoilKnotTube(N, cen)
%PARAMTREFOILKNOTTUBE  A parameterization of the trefoil knot tube.
%   [x,y,z] = paramTrefoilKnotTube(N) returns a mesh with roughly N
%   points in each parameter direction.  SURF(x,y,z) can be used to
%   make a plot.
%
%   [x,y,z] = paramTrefoilKnotTube(N, cen) returns a mesh centered at
%   cen.
%
%   [x,y,z,u,v] = paramTrefoilKnotTube(...) also returns the (u,v)
%   parameter values used (u along the knot, v around the tube's
%   circular cross-section), e.g., for use as the "surfmesh" initial
%   guess in cpParamSurface (see cpTrefoilKnotTube.m).
%
%   The tube is built from a Frenet frame (T,N,B) of the (2,3) torus
%   knot centerline, exactly as in make_trefoil_knot_tube.m (see
%   meshgen/), just evaluated on a plaid (u,v) grid instead of a
%   periodic index grid.

  if (nargin < 1) || isempty(N)
    N = 100;
  end
  if (nargin < 2) || isempty(cen)
    cen = [0 0 0];
  end

  N = max(3, N);

  p = 2;
  q = 3;
  R = 1.65;
  r0 = 0.62;
  rho = 0.19;

  uu = (0:(2*pi/N):(2*pi)).';

  A = R + r0*cos(q*uu);
  Ap = -r0*q*sin(q*uu);
  App = -r0*q^2*cos(q*uu);

  c = [A.*cos(p*uu), A.*sin(p*uu), r0*sin(q*uu)];
  cp = [Ap.*cos(p*uu) - p*A.*sin(p*uu), ...
        Ap.*sin(p*uu) + p*A.*cos(p*uu), ...
        r0*q*cos(q*uu)];
  cpp = [(App - p^2*A).*cos(p*uu) - 2*p*Ap.*sin(p*uu), ...
         (App - p^2*A).*sin(p*uu) + 2*p*Ap.*cos(p*uu), ...
         -r0*q^2*sin(q*uu)];

  T = cp ./ vecnorm(cp, 2, 2);
  nraw = cpp - sum(cpp.*T, 2).*T;
  Nrm = nraw ./ vecnorm(nraw, 2, 2);
  B = cross(T, Nrm, 2);

  vv = 0:(2*pi/N):(2*pi);

  x = c(:,1) + rho*(Nrm(:,1).*cos(vv) + B(:,1).*sin(vv));
  y = c(:,2) + rho*(Nrm(:,2).*cos(vv) + B(:,2).*sin(vv));
  z = c(:,3) + rho*(Nrm(:,3).*cos(vv) + B(:,3).*sin(vv));

  u = repmat(uu, 1, numel(vv));
  v = repmat(vv, numel(uu), 1);

  x = x + cen(1);
  y = y + cen(2);
  z = z + cen(3);
