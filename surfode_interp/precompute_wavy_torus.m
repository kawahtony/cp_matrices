%% Compute closest points of wavy torus

clear ; clc; close all 


%% Load triangular mesh

[vertices, faces] = read_off('data/wavy_torus.off') ;

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
order = 2 ; % Differentiation order

% bw = 1 + sqrt(dim)*(q+1)/2 ; % Bandwidth coefficient
bw = 1.0001*sqrt((dim-1)*((q+1)/2)^2 + ((order/2+(q+1)/2)^2));

xmin = xmin - 0.2 ; xmax = xmax + 0.2 ;
ymin = ymin - 0.2 ; ymax = ymax + 0.2 ;
zmin = zmin - 0.2 ; zmax = zmax + 0.2 ;

x1d = (xmin:dx:xmax)' ; 
y1d = (ymin:dx:ymax)' ;
z1d = (zmin:dx:zmax)' ;

[xx, yy, zz] = meshgrid(x1d, y1d, z1d) ;


%% Compute cp

[cpx, cpy, cpz, dist, uu, vv] = cpWavyTorus(xx, yy, zz) ;


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
    [~, xxu, xxv] = wavy_torus_param(uu_band(k), vv_band(k)) ;
    ntemp = cross(xxu, xxv) ;
    ntemp = ntemp / norm(ntemp) ;
    nx(k) = ntemp(1) ; ny(k) = ntemp(2) ; nz(k) = ntemp(3) ;
end

% figure(1) ; 
% quiver3(cpx_band, cpy_band, cpz_band, nx, ny ,nz)


%% Save

save(['data/wavy_torus_dx=', num2str(dx), '.mat']) 


%% Helper function for computing normal

function [xx, xxu, xxv] = wavy_torus_param(u, v)
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
