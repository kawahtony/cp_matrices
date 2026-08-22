% Driven by mean curvature

clear ; clc ; close all ;


%%

PlyFile = 'bunny.ply';
% PlyFile = 'pig_loop2.ply';
%PlyFile = 'annies_pig.ply';

ptCloud = pcread(PlyFile) ;
Vertices = double(ptCloud.Location) ;
temp = Vertices(:,2) ;
Vertices(:,2) = Vertices(:,3) ; 
Vertices(:,3) = temp ; 
clear temp ;

Faces = MyRobustCrust(Vertices) ;
Faces = double(Faces) ;


%%

% Grid size
dx = 0.025 ;

% make vectors of x, y, positions of the grid
x1d = (-2:dx:2)' ;
y1d = x1d ;
z1d = x1d ;

[xx, yy, zz] = meshgrid(x1d, y1d, z1d) ;


%%

dim = 3 ; % Dimension of embedding space
q = 3 ; % Degree of interpolating polynomial

bw = 1 + sqrt(dim)*(q+1)/2 ;

[sdist, cpx, cpy, cpz, faceidx] = tri2sdf(xx, yy, zz, Faces, Vertices, bw*dx) ;


%% Banding

cpx = cpx(:) ;
cpy = cpy(:) ;
cpz = cpz(:) ;


band = find(abs(sdist) <= bw*dx) ;

% Store closest points in the band
xx_band = xx(band) ;
yy_band = yy(band) ;
zz_band = zz(band) ;
cpx_band = cpx(band) ;
cpy_band = cpy(band) ;
cpz_band = cpz(band) ;
sdist_band = sdist(band) ;


%% Linear operator matrices

Lmat = laplacian_3d_matrix(x1d, y1d, z1d, 2, band, band);
[Dxc, Dyc, Dzc] = firstderiv_cen2_3d_matrices(x1d, y1d, z1d, band, band) ;
Emat = interp3_matrix(x1d, y1d, z1d, cpx_band, cpy_band, cpz_band, q, band);


%% Compute mean curvature and define velocity field

H = 1/2 * Emat*(Lmat*sdist_band) ;
H(H > 10) = 10 ;
H(H < -10) = -10 ;

Hx = Emat * (Dxc * H) ;
Hy = Emat * (Dyc * H) ; 
Hz = Emat * (Dzc * H) ;

%% Plot grids

xp = Vertices(:,1) ;
yp = Vertices(:,2) ;
zp = Vertices(:,3) ;

Eplot = interp3_matrix(x1d, y1d, z1d, xp, yp, zp, 3, band) ;

figure(1) ; clf ; 
trisurf(Faces, xp, yp, zp, Eplot*H, 'FaceAlpha', 0.5) ; 
shading interp
axis equal