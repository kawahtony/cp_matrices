%% Test faster interp

clc ; clear ; close all ;

rng(42)

azure = [0, 128, 255]/255 ;
coral = [255, 127, 80]/255 ;
jade = [0, 168, 107]/255 ;


%% Load closest point of trefoil knot tubular surface

a = 1.5 ; % length along z
b = 1 ; % length along xy plane

dx = 0.05 ;

x1d = (-2.01:dx:2.01)';
y1d = x1d;
z1d = x1d;

[xx, yy, zz] = meshgrid(x1d, y1d, z1d);

[cpx, cpy, cpz, dist] = cpEllipsoid(xx, yy, zz, [a b], [], 'z');


%% Banding

dim = 3 ;
q = 3 ;
bw = 1 + sqrt(dim)*(q + 1)/2 ;

cpx = cpx(:) ; cpy = cpy(:) ; cpz = cpz(:) ;
xx = xx(:) ; yy = yy(:) ; zz = zz(:) ;
dist = dist(:) ;

band = find(dist <= bw*dx);

xx_band = xx(band) ;
yy_band = yy(band) ;
zz_band = zz(band);
cpx_band = cpx(band) ;
cpy_band = cpy(band) ;
cpz_band = cpz(band) ;


phi = pi/3 ;
th = pi/4 ; 
x = [b*sin(phi)*sin(th), b*sin(phi)*cos(th), a*cos(phi)]' ;

tic ;
Etemp = interp3_matrix(x1d, y1d, z1d, x(1), x(2), x(3), 3, band) ;
ans1 = (Etemp*[cpx_band, cpy_band, cpz_band])' ;
toc ;

tic ; 
[Ei, Ej, Es] = interp3_matrix(x1d, y1d, z1d, x(1), x(2), x(3), 3, band) ;
ans2 = [Es'*cpx_band(Ej); Es'*cpy_band(Ej); Es'*cpz_band(Ej)] ;
toc ; 
