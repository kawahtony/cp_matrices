function [H, Hu, Hv] = twistedEllipticTorusMeanCurvature(u,v)
%TWISTEDELLIPTICTORUSMEANCURVATURE  Mean curvature and its u,v partials
%on the twisted elliptic torus.
%   [H, Hu, Hv] = twistedEllipticTorusMeanCurvature(u, v) returns the
%   mean curvature H(u,v) of the three-fold twisted elliptic torus (see
%   cpTwistedEllipticTorus.m for the parametrization X(u,v)), together
%   with the partial derivatives Hu = dH/du and Hv = dH/dv.
%
%   The surface is
%
%     alpha(u)    = 3*u                       (three full twists)
%     radial(u,v) = R + a*cos(v)*cos(alpha) - b*sin(v)*sin(alpha)
%     X(u,v) = ( radial(u,v)*cos(u), radial(u,v)*sin(u), ...
%                a*cos(v)*sin(alpha) + b*sin(v)*cos(alpha) )
%
%   with R=2.15, a=0.62, b=0.34.
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
%   derivatives of X(u,v) (the first two are the same ones computed in
%   cpTwistedEllipticTorus.m's twisted_torus_parm and
%   twisted_torus_parm_2ndpartials helpers, copied verbatim below so
%   this file is self-contained; the third partials are new).
%
%   Verified against direct symbolic differentiation of H(u,v) (sympy)
%   at 50 random (u,v) points; agreement to machine precision.


  sz = size(u);
  u = u(:).';
  v = v(:).';

  [~, Xu, Xv]        = twisted_torus_parm(u, v);
  [Xuu, Xvv, Xuv]     = twisted_torus_parm_2ndpartials(u, v);
  [Xuuu, Xuuv, Xuvv, Xvvv] = twisted_torus_parm_3rdpartials(u, v);

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


end % main twistedEllipticTorusMeanCurvature


%% local helper functions; the first two are copied verbatim from
%% cpTwistedEllipticTorus.m, the third partials were derived the same
%% way (sympy, cross-checked against high-precision finite differences
%% before being transcribed here)

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


function [xxuuu, xxuuv, xxuvv, xxvvv] = twisted_torus_parm_3rdpartials(u, v)
  xxuuu = [ 43*sin(u)/20 + 14*sin(2*u - v)/25 + 48*sin(2*u + v)/25 + 112*sin(4*u - v)/25 + 384*sin(4*u + v)/25;
            -43*cos(u)/20 + 14*cos(2*u - v)/25 + 48*cos(2*u + v)/25 - 112*cos(4*u - v)/25 - 384*cos(4*u + v)/25;
            -189*cos(3*u - v)/50 - 324*cos(3*u + v)/25 ];

  xxuuv = [ -7*sin(2*u - v)/25 + 24*sin(2*u + v)/25 - 28*sin(4*u - v)/25 + 96*sin(4*u + v)/25;
            -7*cos(2*u - v)/25 + 24*cos(2*u + v)/25 + 28*cos(4*u - v)/25 - 96*cos(4*u + v)/25;
            63*cos(3*u - v)/50 - 108*cos(3*u + v)/25 ];

  xxuvv = [ 7*sin(2*u - v)/50 + 12*sin(2*u + v)/25 + 7*sin(4*u - v)/25 + 24*sin(4*u + v)/25;
            7*cos(2*u - v)/50 + 12*cos(2*u + v)/25 - 7*cos(4*u - v)/25 - 24*cos(4*u + v)/25;
            -21*cos(3*u - v)/50 - 36*cos(3*u + v)/25 ];

  xxvvv = [ -7*sin(2*u - v)/100 + 6*sin(2*u + v)/25 - 7*sin(4*u - v)/100 + 6*sin(4*u + v)/25;
            -7*cos(2*u - v)/100 + 6*cos(2*u + v)/25 + 7*cos(4*u - v)/100 - 6*cos(4*u + v)/25;
            7*cos(3*u - v)/50 - 12*cos(3*u + v)/25 ];
end
