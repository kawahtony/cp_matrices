function [cpx,cpy,cpz, dist, uu, vv] = cpTrefoilKnotTube(x, y, z, cen)
%CPTREFOILKNOTTUBE  Closest point function for a tube around a trefoil knot.
%   [cpx,cpy,cpz, dist] = cpTrefoilKnotTube(x,y,z) returns the closest
%   point and distance to (x,y,z) on a smooth circular tube of radius
%   rho=0.19 wrapped around the (2,3) torus (trefoil) knot centerline
%
%     A(u) = R + r0*cos(q*u),   R=1.65, r0=0.62, p=2, q=3
%     c(u) = ( A(u)*cos(p*u), A(u)*sin(p*u), r0*sin(q*u) )
%
%   using the Frenet frame (T,N,B) of c(u):
%
%     X(u,v) = c(u) + rho*( N(u)*cos(v) + B(u)*sin(v) )
%
%   (the same surface as make_trefoil_knot_tube.m, see meshgen/).
%
%   [cpx,cpy,cpz, dist] = cpTrefoilKnotTube(x,y,z, cen) is the same
%   but centered at (xc,yc,zc) = cen.
%
%   [cpx,cpy,cpz, dist, uu, vv] = cpTrefoilKnotTube(...) also returns
%   the (u,v) parameter values of the closest points (u along the
%   knot, v around the tube's cross-section).
%
%   This surface has no closed-form closest point formula, so the
%   closest point is instead found by Newton's method directly on the
%   exact parametric surface (see cpParamSurface.m).
%
%   X(u,v) and its first partials Xu, Xv are exact: because N(u) and
%   B(u) are built from normalized (divided by their own norm)
%   derivatives of c(u), differentiating them by hand is impractical,
%   so Xu, Xv below were derived with symbolic differentiation (sympy)
%   and cross-checked against high-precision (50-digit) finite
%   differences before being transcribed here -- hence the terse
%   s0, s1, ... intermediate variables (common subexpressions shared
%   by X, Xu and Xv).  The Newton solve only needs the *second*
%   partials Xuu, Xvv, Xuv as a curvature estimate for its step, not
%   for accuracy (the reported closest point/distance only ever use X
%   and the converged (u,v)), so rather than transcribe the far larger
%   exact closed form for the second partials, they are computed
%   below with a central finite difference of the (already exact) Xu,
%   Xv -- simple, robust, and easy to check by eye.
%
%   Note: unlike most other cp functions here, this does not return a
%   signed distance (there's no cheap/robust inside test for a knotted
%   tube); dist is always >= 0.

  % defaults
  if (nargin < 4) || isempty(cen)
    cen = [0 0 0];
  end

  % shift to the origin
  x = x - cen(1);
  y = y - cen(2);
  z = z - cen(3);

  paramf = @trefoil_parm;
  paramf2nd = @trefoil_parm_2ndpartials;
  paramfEdge = [];  % no boundary: this is a closed surface
  optf = [];
  surfmesh = {[] [] [] [] []};
  [surfmesh{:}] = paramTrefoilKnotTube(160);
  paramAdjust = [];

  How = 2;  % Newton's method (fastest, see cpParamSurface.m)
  [cpx,cpy,cpz,dist,~,uu,vv] = cpParamSurface(x,y,z, paramf, paramf2nd, paramfEdge, optf, surfmesh, [-inf -inf], [inf inf], paramAdjust, How);

  % shift back
  cpx = cpx + cen(1);
  cpy = cpy + cen(2);
  cpz = cpz + cen(3);

end % main cpTrefoilKnotTube


%% local helper functions

function [xx, xxu, xxv] = trefoil_parm(u, v)
%TREFOIL_PARM  Position and exact first partials of the tube X(u,v).
%   Derived with sympy (symbolic differentiation of the Frenet-frame
%   construction used by make_trefoil_knot_tube.m) and verified
%   against 50-digit finite differences; see cpTrefoilKnotTube.m.

  s0 = 2*u;
  s1 = cos(s0);
  s2 = 3*u;
  s3 = cos(s2);
  s4 = 62*s3 + 165;
  s5 = s1.*s4;
  s6 = 2*s5;
  s7 = sin(s0);
  s8 = sin(s2);
  s9 = s7.*s8;
  s10 = s1.*s3;
  s11 = 4*s5;
  s12 = s1.*s8;
  s13 = s4.*s7;
  s14 = 93*s12 + s13;
  s15 = cos(u);
  s16 = s15.^2;
  s17 = s15.^3;
  s18 = 30752*s15.^6 - 46128*s15.^4 - 30690*s15 + 17298*s16 + 40920*s17 + 17937;
  s19 = 1./s18;
  s20 = 372*s9;
  s21 = 279*s10 - s20 + s6;
  s22 = -s21;
  s23 = s5 - 93*s9;
  s24 = s3.*s7;
  s25 = 372*s12;
  s26 = 2*s13;
  s27 = 279*s24 + s25 + s26;
  s28 = s14.*s22 + s23.*s27 + 25947*s3.*s8;
  s29 = s19.*s28;
  s30 = 558*s10 + s11 + s14.*s29 - 744*s9;
  s31 = s8.^2;
  s32 = sin(u);
  s33 = 10230*s3 - 30752*s32.^6 + 46128*s32.^4 - 17298*s32.^2 + 19859;
  s34 = 1./s33;
  s35 = s3.*s34;
  s36 = 31*s4;
  s37 = s35.*s36 - 1;
  s38 = 31*s32;
  s39 = 5*u;
  s40 = sin(s39);
  s41 = s38 + 155*s40 + 330*s7;
  s42 = 165*s8;
  s43 = 6*u;
  s44 = s42 + 31*sin(s43);
  s45 = 93*s44;
  s46 = s34.*s45;
  s47 = 31*s15;
  s48 = cos(s39);
  s49 = 660*s1 + s47 + 775*s48;
  s50 = s41.*s46 + s49;
  s51 = 660*s7;
  s52 = 775*s40;
  s53 = -s47;
  s54 = 330*s1 + 155*s48 + s53;
  s55 = s34.*s54;
  s56 = s38 + s45.*s55 - s51 - s52;
  s57 = 311364*s31.*s37.^2 + s50.^2 + s56.^2;
  s58 = 1./sqrt(s57);
  s59 = cos(v);
  s60 = s58.*s59;
  s61 = 38*s60;
  s62 = s29.*s3 - 6*s8;
  s63 = -s62;
  s64 = s23.*s63;
  s65 = 4*s13;
  s66 = 744*s12 - s23.*s29 + 558*s24 + s65;
  s67 = s3.*s66;
  s68 = s64 - s67;
  s69 = sin(v);
  s70 = sqrt(2);
  s71 = 1./sqrt(s18);
  s72 = s70.*s71;
  s73 = s58.*s69.*s72;
  s74 = 1767*s73;
  s75 = s14.*s62;
  s76 = -s30;
  s77 = s3.*s76;
  s78 = s75 + s77;
  s79 = 3534*s60;
  s80 = s23.*s30;
  s81 = s14.*s66;
  s82 = s80 + s81;
  s83 = 19*s73;
  s84 = s3.^2;
  s85 = s14.*(1953*s12 + 1674*s24 + s65) + s23.*(1674*s10 + s11 - 1953*s9) - s27.^2 - 77841*s31 + 77841*s84;
  s86 = s19.*(s22.^2 - s85);
  s87 = 186*s15;
  s88 = 992*s15.^5 + 660*s16 - 992*s17 + s87 - 165;
  s89 = s18.^(-2);
  s90 = 3906*s12 + 8*s13 - 186*s14.*s28.*s32.*s88.*s89 + 3348*s24;
  s91 = s14.*s86 + s22.*s29 + s90;
  s92 = s57.^(-3/2);
  s93 = s33.^(-2);
  s94 = 992*s15;
  s95 = s32.^5.*s94 - s32.^3.*s94 + s32.*s87 + s42;
  s96 = -1922*s3.*s4.*s8.*s93.*s95 + s3 + s31.*s34.*s36 + 1922*s31.*s35 - 31*s34.*s4.*s84;
  s97 = 934092*s37.*s8;
  s98 = 1320*s1;
  s99 = 3875*s48;
  s100 = -31*s32;
  s101 = s100 + s51 + s52;
  s102 = 165*s3 + 62*cos(s43);
  s103 = 279*s102;
  s104 = 17298*s44.*s93.*s95;
  s105 = s100 + s103.*s34.*s41 + s104.*s41 - 3875*s40 + s46.*s49 - 1320*s7;
  s106 = s105.*s50 + s56.*(-s101.*s46 + s103.*s55 + s104.*s54 + s47 - s98 - s99) - s96.*s97;
  s107 = 3*s8;
  s108 = s32.*s88;
  s109 = 93*s108.*s19;
  s110 = 186*s108.*s28.*s89;
  s111 = s107.*s29 - s110.*s3 + 18*s3;
  s112 = s111 + s3.*s86;
  s113 = 3348*s10 - s110.*s23 + s27.*s29 + 8*s5 - 3906*s9;
  s114 = s113 + s23.*s86;
  s115 = 1./s57;
  s116 = -s106.*s115;
  s117 = s19.*(s21.^2 - s85);
  s118 = s105.*s50 + s56.*(-s101.*s46 + 279*s102.*s34.*s54 + 17298*s44.*s54.*s93.*s95 - s53 - s98 - s99) - s96.*s97;
  s119 = s115.*s118;
  s120 = s59.*s72;
  s121 = 93*s120;
  s122 = 19*s58/200;

  xx  = [ -s30.*s61/200 + s6/200 - s68.*s74/200;
          s26/200 - s61.*s66/200 + s74.*s78/200;
          -s63.*s79/200 + 31*s8/50 + s82.*s83/200 ];

  xxu = [ 19*s106.*s30.*s59.*s92/100 - s25/200 + 1767*s58.*s69.*s70.*s71.*(s106.*s115.*s23.*s63 - s106.*s115.*s67 - s107.*s66 - s109.*s64 - s112.*s23 + s114.*s3 + 93*s19.*s3.*s32.*s66.*s88 + s27.*s63)/200 + s61.*s91/200 - s65/200;
          19*s106.*s59.*s66.*s92/100 + s11/200 - s114.*s61/200 - s20/200 - s74.*(-s109.*s75 - s109.*s77 + s112.*s14 - s116.*s75 - s116.*s77 + s22.*s62 - s3.*s91 + 3*s76.*s8)/200;
          1767*s118.*s59.*s63.*s92/100 + 93*s3/50 - s79.*(s111 + s117.*s3)/200 + s83.*(s109.*s80 + s109.*s81 - s119.*s80 - s119.*s81 + s14.*(s113 + s117.*s23) + s21.*s66 + s23.*(-s117.*s14 + s19.*s21.*s28 - s90) - s27.*s30)/200 ];

  xxv = [ s122.*(-s121.*s68 + 2*s30.*s69);
          s122.*(s121.*s78 + 2*s66.*s69);
          s122.*(s120.*s82 + 186*s63.*s69) ];
end


function [xxuu, xxvv, xxuv] = trefoil_parm_2ndpartials(u, v)
%TREFOIL_PARM_2NDPARTIALS  Approximate second partials of X(u,v).
%   Xu and Xv from trefoil_parm are exact; Newton's method (see
%   cpParamSurface.m / helper_newton) only uses the second partials to
%   build a local quadratic model for its step direction, not for
%   accuracy of the final answer, so a plain central finite difference
%   of the exact first partials is used here instead of the (much
%   larger and harder to transcribe correctly) exact closed form.

  h = 1e-6;

  [~, xxu_up, ~]         = trefoil_parm(u+h, v);
  [~, xxu_um, ~]         = trefoil_parm(u-h, v);
  [~, xxu_vp, xxv_vp]    = trefoil_parm(u, v+h);
  [~, xxu_vm, xxv_vm]    = trefoil_parm(u, v-h);

  xxuu = (xxu_up - xxu_um) / (2*h);
  xxvv = (xxv_vp - xxv_vm) / (2*h);
  xxuv = (xxu_vp - xxu_vm) / (2*h);
end
