%% Load sdf of frog surface directly

clear ; clc ; close all

jade = [0, 168, 107]/255 ;

load('data/frog_dx=0.025.mat')

[face, vertex] = isosurface(xx, yy, zz, sdist, 0) ;

sdist = sdist(:) ;
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

xp = vertex(:,1) ;
yp = vertex(:,2) ;
zp = vertex(:,3) ;

Eplot = interp3_matrix(x1d, y1d, z1d, xp, yp, zp, 3, band) ;

figure(1) ; clf ; 
trisurf(face, xp, yp, zp, Eplot*H, 'FaceAlpha', 0.5) ; 
shading interp
axis equal
