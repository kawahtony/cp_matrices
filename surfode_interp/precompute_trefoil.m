%% Compute closest point function of trefoil knot tube


clear ; clc; close all 


%% Load triangular mesh

surface = 'trefoil_knot_tube' ;

[vertices, faces] = read_off(['data/', surface, '.off']) ;

figure(1) ; clf ; hold on ;
trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3), 'FaceColor',[0.2 0.5 0.9], 'EdgeColor','none');
axis equal
camlight
lighting gouraud

temp = xlim ; xmin = temp(1) ; xmax = temp(2) ;
temp = ylim ; ymin = temp(1) ; ymax = temp(2) ;
temp = zlim ; zmin = temp(1) ; zmax = temp(2) ;


%% Cartesian mesh in the embedding space

dx = 0.025 ; % Grid size

dim = 3 ;  % Dimension of embedding space
q = 3 ; % Degree of interpolating polynomial

bw = 1 + sqrt(dim)*(q+1)/2 ; % Bandwidth coefficient

xmin = xmin - 0.2 ; xmax = xmax + 0.2 ;
ymin = ymin - 0.2 ; ymax = ymax + 0.2 ;
zmin = zmin - 0.2 ; zmax = zmax + 0.2 ;

x1d = (xmin:dx:xmax)' ; 
y1d = (ymin:dx:ymax)' ;
z1d = (zmin:dx:zmax)' ;

[xx, yy, zz] = meshgrid(x1d, y1d, z1d) ;


%% Compute cp

[cpx, cpy, cpz, dist, uu, vv] = cpTrefoilKnotTube(xx, yy, zz) ;


%% Banding

xg = xx(:) ;
yg = yy(:) ;
zg = zz(:) ;
cpxg = cpx(:) ;
cpyg = cpy(:) ;
cpzg = cpz(:) ;
distg = dist(:) ;
uug = uu(:) ;
vvg = vv(:) ;

band = find(abs(dist) <= bw*dx) ;

% Store closest points in the band
xx_band = xg(band) ;
yy_band = yg(band) ;
zz_band = zg(band) ;
cpx_band = cpxg(band) ;
cpy_band = cpyg(band) ;
cpz_band = cpzg(band) ;
dist_band = distg(band) ;
uu_band = uug(band) ;
vv_band = vvg(band) ;

nx = zeros(size(band)) ; ny = nx ; nz = nx ;

for k = 1:length(band)
    [~, xxu, xxv] = trefoil_parm(uu_band(k), vv_band(k)) ;
    ntemp = cross(xxu, xxv) ;
    ntemp = ntemp / norm(ntemp) ;
    nx(k) = ntemp(1) ; ny(k) = ntemp(2) ; nz(k) = ntemp(3) ;
end


% figure(1) ; 
% quiver3(cpx_band, cpy_band, cpz_band, nx, ny ,nz)


%% Save

save(['data/', surface, '_dx=', num2str(dx), '.mat']) 

%%

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