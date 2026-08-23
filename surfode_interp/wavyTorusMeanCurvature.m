function [H, Hu, Hv] = wavyTorusMeanCurvature(u,v)
%WAVYTORUSMEANCURVATURE  Mean curvature and its u,v partials on the
%wavy torus.
%   [H, Hu, Hv] = wavyTorusMeanCurvature(u, v) returns the mean
%   curvature H(u,v) of the five-fold corrugated torus (see
%   cpWavyTorus.m for the parametrization X(u,v)), together with the
%   partial derivatives Hu = dH/du and Hv = dH/dv.
%
%   u, v may be scalars or arrays of the same size; H, Hu, Hv are
%   returned in the same shape.
%
%   Method: standard first/second fundamental form formula
%
%     H = (E*N - 2*F*M + G*L) / (2*(E*G - F^2))
%
%   with E=Xu.Xu, F=Xu.Xv, G=Xv.Xv, unit normal n=(Xu x Xv)/|Xu x Xv|,
%   and L=Xuu.n, M=Xuv.n, N=Xvv.n.  Hu, Hv are obtained by
%   differentiating this formula (product/quotient rule), which in
%   turn requires the derivative of the unit normal n; this is done
%   analytically using the exact first, second, and third partial
%   derivatives of X(u,v) (the same ones computed in cpWavyTorus.m's
%   wavy_torus_parm, wavy_torus_parm_2ndpartials, and
%   wavy_torus_parm_3rdpartials helpers, copied verbatim below so this
%   file is self-contained).
%
%   Verified against direct symbolic differentiation of H(u,v) (sympy)
%   at several test points; agreement to machine precision.


  sz = size(u);
  u = u(:).';
  v = v(:).';

  [~, Xu, Xv]        = wavy_torus_parm(u, v);
  [Xuu, Xvv, Xuv]     = wavy_torus_parm_2ndpartials(u, v);
  [Xuuu, Xuuv, Xuvv, Xvvv] = wavy_torus_parm_3rdpartials(u, v);

  % first fundamental form
  E = dot(Xu, Xu);
  F = dot(Xu, Xv);
  G = dot(Xv, Xv);

  % unit normal
  Ncross = cross(Xu, Xv);
  Nnorm  = sqrt(dot(Ncross, Ncross));
  n      = Ncross ./ Nnorm;

  % second fundamental form
  L = dot(Xuu, n);
  M = dot(Xuv, n);
  N = dot(Xvv, n);

  num = E.*N - 2*F.*M + G.*L;
  den = 2*(E.*G - F.^2);
  H   = num ./ den;

  % derivatives of E, F, G
  Eu = 2*dot(Xu, Xuu);
  Ev = 2*dot(Xu, Xuv);
  Fu = dot(Xuu, Xv) + dot(Xu, Xuv);
  Fv = dot(Xuv, Xv) + dot(Xu, Xvv);
  Gu = 2*dot(Xv, Xuv);
  Gv = 2*dot(Xv, Xvv);

  % derivative of the (unnormalized) cross-product normal
  dNcross_du = cross(Xuu, Xv) + cross(Xu, Xuv);
  dNcross_dv = cross(Xuv, Xv) + cross(Xu, Xvv);

  % derivative of the unit normal (tangential projection, divided by norm)
  nu = (dNcross_du - n.*dot(n, dNcross_du)) ./ Nnorm;
  nv = (dNcross_dv - n.*dot(n, dNcross_dv)) ./ Nnorm;

  % derivatives of L, M, N
  Lu = dot(Xuuu, n) + dot(Xuu, nu);
  Lv = dot(Xuuv, n) + dot(Xuu, nv);
  Mu = dot(Xuuv, n) + dot(Xuv, nu);
  Mv = dot(Xuvv, n) + dot(Xuv, nv);
  Nu = dot(Xuvv, n) + dot(Xvv, nu);
  Nv = dot(Xvvv, n) + dot(Xvv, nv);

  % quotient rule on H = num/den
  num_u = Eu.*N + E.*Nu - 2*(Fu.*M + F.*Mu) + Gu.*L + G.*Lu;
  num_v = Ev.*N + E.*Nv - 2*(Fv.*M + F.*Mv) + Gv.*L + G.*Lv;

  den_u = 2*(Eu.*G + E.*Gu - 2*F.*Fu);
  den_v = 2*(Ev.*G + E.*Gv - 2*F.*Fv);

  Hu = (num_u.*den - num.*den_u) ./ den.^2;
  Hv = (num_v.*den - num.*den_v) ./ den.^2;

  H  = reshape(H,  sz);
  Hu = reshape(Hu, sz);
  Hv = reshape(Hv, sz);


end % main wavyTorusMeanCurvature


%% local helper functions, copied verbatim from cpWavyTorus.m (surface
%% derived and cross-checked with sympy; formulas verified against
%% high-precision finite differences before being transcribed here)

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


function [xxuuu, xxuuv, xxuvv, xxvvv] = wavy_torus_parm_3rdpartials(u, v)
  xxuuu = [ 41*sin(u)/20 + 192*sin(4*u)/25 + 648*sin(6*u)/25 + 13*sin(u - v)/50 + 13*sin(u + v)/50 + 169*sin(u + 2*v)/10000 + 169*sin(u + 4*v)/10000 + 4563*sin(3*u + 2*v)/10000 + 4563*sin(3*u + 4*v)/10000;
            -41*cos(u)/20 + 192*cos(4*u)/25 - 648*cos(6*u)/25 - 13*cos(u - v)/50 - 13*cos(u + v)/50 + 169*cos(u + 2*v)/10000 + 169*cos(u + 4*v)/10000 - 4563*cos(3*u + 2*v)/10000 - 4563*cos(3*u + 4*v)/10000;
            338*sin(v).*sin(2*u + 3*v)/625 - 45*cos(5*u)/2 ];

  xxuuv = [ -13*sin(u - v)/50 + 13*sin(u + v)/50 + 169*sin(u + 2*v)/5000 + 169*sin(u + 4*v)/2500 + 1521*sin(3*u + 2*v)/5000 + 1521*sin(3*u + 4*v)/2500;
            13*cos(u - v)/50 - 13*cos(u + v)/50 + 169*cos(u + 2*v)/5000 + 169*cos(u + 4*v)/2500 - 1521*cos(3*u + 2*v)/5000 - 1521*cos(3*u + 4*v)/2500;
            169*cos(2*u + 2*v)/625 - 338*cos(2*u + 4*v)/625 ];

  xxuvv = [ 13*sin(u - v)/50 + 13*sin(u + v)/50 + 169*sin(u + 2*v)/2500 + 169*sin(u + 4*v)/625 + 507*sin(3*u + 2*v)/2500 + 507*sin(3*u + 4*v)/625;
            -13*cos(u - v)/50 - 13*cos(u + v)/50 + 169*cos(u + 2*v)/2500 + 169*cos(u + 4*v)/625 - 507*cos(3*u + 2*v)/2500 - 507*cos(3*u + 4*v)/625;
            169*cos(2*u + 2*v)/625 - 676*cos(2*u + 4*v)/625 ];

  xxvvv = [ 13*(25*sin(v) + 13*sin(2*u + 2*v) + 104*sin(2*u + 4*v)).*cos(u)/625;
            13*(25*sin(v) + 13*sin(2*u + 2*v) + 104*sin(2*u + 4*v)).*sin(u)/625;
            -13*cos(v)/25 + 169*cos(2*u + 2*v)/625 - 1352*cos(2*u + 4*v)/625 ];
end