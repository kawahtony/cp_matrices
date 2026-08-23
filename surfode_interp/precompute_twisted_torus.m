%% Compute closest point of twisted elliptic torus

clear ; clc; close all 


%% Cartesian mesh in the embedding space

dx = 0.025 ; % Grid size

dim = 3 ;  % Dimension of embedding space
q = 3 ; % Degree of interpolating polynomial
order = 2 ; % Differentiation order

% bw = 1 + sqrt(dim)*(q+1)/2 ; % Bandwidth coefficient
bw = 1.0001*sqrt((dim-1)*((q+1)/2)^2 + ((order/2+(q+1)/2)^2));

xmin = -4.940046325278577 ; xmax = 4.940046325278577 ;
ymin = -0.32 ; ymax = 0.32 ;
zmin = -0.82 ; zmax = 0.82 ;

x1d = (xmin:dx:xmax)' ; 
y1d = (ymin:dx:ymax)' ;
z1d = (zmin:dx:zmax)' ;

[xx, yy, zz] = meshgrid(x1d, y1d, z1d) ;


%% Compute cp

[cpx, cpy, cpz, dist, uu, vv] = cpTwistedEllipticTorus(xx, yy, zz) ;


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


%% Compute unit outward normal

nx = zeros(size(band)) ; ny = nx ; nz = nx ;

for k = 1:length(band)
    [~, xxu, xxv] = twisted_torus_param(uu_band(k), vv_band(k)) ;
    ntemp = cross(xxu, xxv) ;
    ntemp = ntemp / norm(ntemp) ;
    nx(k) = ntemp(1) ; ny(k) = ntemp(2) ; nz(k) = ntemp(3) ;
end


% figure(1) ; 
% quiver3(cpx_band, cpy_band, cpz_band, nx, ny ,nz)


%% Save

save(['data/twisted_elliptic_torus_dx=', num2str(dx), '.mat']) 


%% Helpful function for computing normal


function [xx, xxu, xxv] = twisted_torus_param(u, v)
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
