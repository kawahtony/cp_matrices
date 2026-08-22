function [cpx,cpy,cpz, dist, uu, vv] = cpWavyTorus(x, y, z, cen)
%CPWAVYTORUS  Closest point function for a five-fold corrugated torus.
%   [cpx,cpy,cpz, dist] = cpWavyTorus(x,y,z) returns the closest point
%   and distance to (x,y,z) on the wavy torus
%
%     R(u)  = 2.05 + 0.24*cos(5u)
%     zc(u) = 0.18*sin(5u)
%     r(u,v) = 0.52*(1 + 0.13*cos(3v+2u))
%     X(u,v) = ( (R(u)+r(u,v)*cos(v))*cos(u), ...
%                (R(u)+r(u,v)*cos(v))*sin(u), ...
%                zc(u) + r(u,v)*sin(v) )
%
%   (the same surface as make_wavy_torus.m, see meshgen/).
%
%   [cpx,cpy,cpz, dist] = cpWavyTorus(x,y,z, cen) is the same but
%   centered at (xc,yc,zc) = cen.
%
%   [cpx,cpy,cpz, dist, uu, vv] = cpWavyTorus(...) also returns the
%   (u,v) parameter values of the closest points.
%
%   Unlike a plain torus this surface has no simple closed-form
%   closest point formula (R, zc, r are all functions of u, and r
%   also depends on v), so the closest point is instead found by
%   Newton's method directly on the exact parametric surface (see
%   cpParamSurface.m), using the analytic first and second partial
%   derivatives of X(u,v) below.  This is exact (up to the Newton
%   solve's tolerance), unlike a mesh/triangulation-based closest
%   point.
%
%   Note: unlike most other cp functions here, this does not return a
%   signed distance (there's no cheap/robust inside test for this
%   corrugated shape); dist is always >= 0.

  % defaults
  if (nargin < 4) || isempty(cen)
    cen = [0 0 0];
  end

  % shift to the origin
  x = x - cen(1);
  y = y - cen(2);
  z = z - cen(3);

  paramf = @wavy_torus_parm;
  paramf2nd = @wavy_torus_parm_2ndpartials;
  paramfEdge = [];  % no boundary: this is a closed surface
  optf = [];
  surfmesh = {[] [] [] [] []};
  [surfmesh{:}] = paramWavyTorus(120);
  paramAdjust = [];

  How = 2;  % Newton's method (fastest, see cpParamSurface.m)
  [cpx,cpy,cpz,dist,~,uu,vv] = cpParamSurface(x,y,z, paramf, paramf2nd, paramfEdge, optf, surfmesh, [-inf -inf], [inf inf], paramAdjust, How);

  % shift back
  cpx = cpx + cen(1);
  cpy = cpy + cen(2);
  cpz = cpz + cen(3);

end % main cpWavyTorus


%% local helper functions (surface derived and cross-checked with sympy,
%% see cp_matrices dev notes; formulas verified against high-precision
%% finite differences before being transcribed here)

function [xx, xxu, xxv] = wavy_torus_parm(u, v)
  xx  = [ ((169*cos(2*u + 3*v)/2500 + 13/25).*cos(v) + 6*cos(5*u)/25 + 41/20).*cos(u);
          ((169*cos(2*u + 3*v)/2500 + 13/25).*cos(v) + 6*cos(5*u)/25 + 41/20).*sin(u);
          (169*cos(2*u + 3*v)/2500 + 13/25).*sin(v) + 9*sin(5*u)/50 ];

  xxu = [ -(1500*sin(5*u) + 169*sin(2*u + 3*v).*cos(v)).*cos(u)/1250 - (13*(13*cos(2*u + 3*v) + 100).*cos(v) + 600*cos(5*u) + 5125).*sin(u)/2500;
          -(1500*sin(5*u) + 169*sin(2*u + 3*v).*cos(v)).*sin(u)/1250 + (13*(13*cos(2*u + 3*v) + 100).*cos(v) + 600*cos(5*u) + 5125).*cos(u)/2500;
          -169*sin(v).*sin(2*u + 3*v)/1250 + 9*cos(5*u)/10 ];

  xxv = [ -13*(100*sin(v) + 13*sin(2*u + 2*v) + 26*sin(2*u + 4*v)).*cos(u)/2500;
          -13*(100*sin(v) + 13*sin(2*u + 2*v) + 26*sin(2*u + 4*v)).*sin(u)/2500;
          13*cos(v)/25 - 169*cos(2*u + 2*v)/2500 + 169*cos(2*u + 4*v)/1250 ];
end


function [xxuu, xxvv, xxuv] = wavy_torus_parm_2ndpartials(u, v)
  xxuu = [ -41*cos(u)/20 - 48*cos(4*u)/25 - 108*cos(6*u)/25 - 13*cos(u - v)/50 - 13*cos(u + v)/50 - 169*cos(u + 2*v)/10000 - 169*cos(u + 4*v)/10000 - 1521*cos(3*u + 2*v)/10000 - 1521*cos(3*u + 4*v)/10000;
           -41*sin(u)/20 + 48*sin(4*u)/25 - 108*sin(6*u)/25 - 13*sin(u - v)/50 - 13*sin(u + v)/50 + 169*sin(u + 2*v)/10000 + 169*sin(u + 4*v)/10000 - 1521*sin(3*u + 2*v)/10000 - 1521*sin(3*u + 4*v)/10000;
           -9*sin(5*u)/2 - 169*sin(v).*cos(2*u + 3*v)/625 ];

  xxvv = [ -13*(50*cos(v) + 13*cos(2*u + 2*v) + 52*cos(2*u + 4*v)).*cos(u)/1250;
           -13*(50*cos(v) + 13*cos(2*u + 2*v) + 52*cos(2*u + 4*v)).*sin(u)/1250;
           -13*sin(v)/25 + 169*sin(2*u + 2*v)/1250 - 338*sin(2*u + 4*v)/625 ];

  xxuv = [ -169*(cos(2*u + 2*v) + 2*cos(2*u + 4*v)).*cos(u)/1250 + 13*(100*sin(v) + 13*sin(2*u + 2*v) + 26*sin(2*u + 4*v)).*sin(u)/2500;
           -169*(cos(2*u + 2*v) + 2*cos(2*u + 4*v)).*sin(u)/1250 - 13*(100*sin(v) + 13*sin(2*u + 2*v) + 26*sin(2*u + 4*v)).*cos(u)/2500;
           169*sin(2*u + 2*v)/1250 - 169*sin(2*u + 4*v)/625 ];
end
