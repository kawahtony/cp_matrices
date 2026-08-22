function [cpx,cpy,cpz, dist, uu, vv] = cpTwistedEllipticTorus(x, y, z, cen)
%CPTWISTEDELLIPTICTORUS  Closest point function for a twisted elliptic torus.
%   [cpx,cpy,cpz, dist] = cpTwistedEllipticTorus(x,y,z) returns the
%   closest point and distance to (x,y,z) on the surface
%
%     alpha(u)    = 3*u                       (three full twists)
%     radial(u,v) = R + a*cos(v)*cos(alpha) - b*sin(v)*sin(alpha)
%     X(u,v) = ( radial(u,v)*cos(u), radial(u,v)*sin(u), ...
%                a*cos(v)*sin(alpha) + b*sin(v)*cos(alpha) )
%
%   with R=2.15, a=0.62, b=0.34 (the same surface as
%   make_twisted_elliptic_torus.m, see meshgen/).
%
%   [cpx,cpy,cpz, dist] = cpTwistedEllipticTorus(x,y,z, cen) is the
%   same but centered at (xc,yc,zc) = cen.
%
%   [cpx,cpy,cpz, dist, uu, vv] = cpTwistedEllipticTorus(...) also
%   returns the (u,v) parameter values of the closest points.
%
%   Because the tube's elliptical cross-section rotates as it travels
%   around the loop, this surface has no simple closed-form closest
%   point formula, so the closest point is instead found by Newton's
%   method directly on the exact parametric surface (see
%   cpParamSurface.m), using the analytic first and second partial
%   derivatives of X(u,v) below.  This is exact (up to the Newton
%   solve's tolerance), unlike a mesh/triangulation-based closest
%   point.
%
%   Note: unlike most other cp functions here, this does not return a
%   signed distance (there's no cheap/robust inside test for this
%   twisted shape); dist is always >= 0.

  % defaults
  if (nargin < 4) || isempty(cen)
    cen = [0 0 0];
  end

  % shift to the origin
  x = x - cen(1);
  y = y - cen(2);
  z = z - cen(3);

  paramf = @twisted_torus_parm;
  paramf2nd = @twisted_torus_parm_2ndpartials;
  paramfEdge = [];  % no boundary: this is a closed surface
  optf = [];
  surfmesh = {[] [] [] [] []};
  [surfmesh{:}] = paramTwistedEllipticTorus(120);
  paramAdjust = [];

  How = 2;  % Newton's method (fastest, see cpParamSurface.m)
  [cpx,cpy,cpz,dist,~,uu,vv] = cpParamSurface(x,y,z, paramf, paramf2nd, paramfEdge, optf, surfmesh, [-inf -inf], [inf inf], paramAdjust, How);

  % shift back
  cpx = cpx + cen(1);
  cpy = cpy + cen(2);
  cpz = cpz + cen(3);

end % main cpTwistedEllipticTorus


%% local helper functions (derived and cross-checked with sympy; formulas
%% verified against high-precision finite differences before being
%% transcribed here)

function [xx, xxu, xxv] = twisted_torus_parm(u, v)
  xx  = [ (-17*sin(3*u).*sin(v)/50 + 31*cos(3*u).*cos(v)/50 + 43/20).*cos(u);
          (-17*sin(3*u).*sin(v)/50 + 31*cos(3*u).*cos(v)/50 + 43/20).*sin(u);
          31*sin(3*u).*cos(v)/50 + 17*sin(v).*cos(3*u)/50 ];

  xxu = [ -43*sin(u)/20 - 7*sin(2*u - v)/50 - 12*sin(2*u + v)/25 - 7*sin(4*u - v)/25 - 24*sin(4*u + v)/25;
          43*cos(u)/20 - 7*cos(2*u - v)/50 - 12*cos(2*u + v)/25 + 7*cos(4*u - v)/25 + 24*cos(4*u + v)/25;
          21*cos(3*u - v)/50 + 36*cos(3*u + v)/25 ];

  xxv = [ (7*sin(3*u - v) - 24*sin(3*u + v)).*cos(u)/50;
          (7*sin(3*u - v) - 24*sin(3*u + v)).*sin(u)/50;
          -7*cos(3*u - v)/50 + 12*cos(3*u + v)/25 ];
end


function [xxuu, xxvv, xxuv] = twisted_torus_parm_2ndpartials(u, v)
  xxuu = [ -43*cos(u)/20 - 7*cos(2*u - v)/25 - 24*cos(2*u + v)/25 - 28*cos(4*u - v)/25 - 96*cos(4*u + v)/25;
           -43*sin(u)/20 + 7*sin(2*u - v)/25 + 24*sin(2*u + v)/25 - 28*sin(4*u - v)/25 - 96*sin(4*u + v)/25;
           -63*sin(3*u - v)/50 - 108*sin(3*u + v)/25 ];

  xxvv = [ -(7*cos(3*u - v) + 24*cos(3*u + v)).*cos(u)/50;
           -(7*cos(3*u - v) + 24*cos(3*u + v)).*sin(u)/50;
           -7*sin(3*u - v)/50 - 12*sin(3*u + v)/25 ];

  xxuv = [ 7*cos(2*u - v)/50 - 12*cos(2*u + v)/25 + 7*cos(4*u - v)/25 - 24*cos(4*u + v)/25;
           -7*sin(2*u - v)/50 + 12*sin(2*u + v)/25 + 7*sin(4*u - v)/25 - 24*sin(4*u + v)/25;
           21*sin(3*u - v)/50 - 36*sin(3*u + v)/25 ];
end
