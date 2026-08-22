% Test: compute mean curvature

clear ; clc ; close all ;


%% Discretizaiton of embedding space

% Mesh size
dx = 0.05 ;

% Interpolation order
q = 3 ; 

% Make vectors of x, y, positions of the grid
x1d = (-2:dx:2)' ;
y1d = x1d ;
z1d = x1d ;

% Embedding space
[xx, yy, zz] = meshgrid(x1d, y1d, z1d) ;

% Compute closest point 
R = 1 ;
r = 0.4 ;
[cpx, cpy, cpz, sdist] = cpTorus(xx, yy, zz, R, r) ;
[faces, vertices] = isosurface(xx, yy, zz, sdist, 0) ;
xx = xx(:) ; 
yy = yy(:) ; 
zz = zz(:) ;
cpx = cpx(:) ;
cpy = cpy(:) ;
cpz = cpz(:) ;
sdist = sdist(:) ;

% Computational band
dim = 3 ;
% bw = 2*sqrt(dim) ;
bw = sqrt(dim*((q + 1)/2)^2) ;
band = find(abs(sdist) <= bw*dx) ;

% Store relevant closest point information in the band
xx_band = xx(band) ;
yy_band = yy(band) ; 
zz_band = zz(band) ;
cpx_band = cpx(band) ;
cpy_band = cpy(band) ;
cpz_band = cpz(band) ;
sdist_band = sdist(band) ;



% phi_band = atan2(cpy_band, cpx_band) ;
% th_band = acos(-cpz_band/r) ;
% 
% H_exact = 1/2 * (1/r + sin(th_band)./(R + r*sin(th_band)) ) ;


%% Construct an interpolation matrix for plotting on surface

% plotting grid on torus, based on parametrization
[xp,yp,zp] = paramTorus(200, R, r) ;
xp1 = xp(:) ; yp1 = yp(:) ; zp1 = zp(:) ;

% % Eplot is a matrix which interpolates data onto the plotting grid
Eplot = interp3_matrix(x1d, y1d, z1d, xp1, yp1, zp1, 3, band) ;



%% Linear operator matrices


Lmat = laplacian_3d_matrix(x1d, y1d, z1d, 2, band, band);
[Dxc, Dyc, Dzc] = firstderiv_cen2_3d_matrices(x1d, y1d, z1d, band, band) ;
Emat = interp3_matrix(x1d, y1d, z1d, cpx_band, cpy_band, cpz_band, q, band);


%% Compute mean curvature and define velocity field

H = 1/2 * Emat*(Lmat*sdist_band) ;

Hx = Emat * (Dxc * H) ;
Hy = Emat * (Dyc * H) ; 
Hz = Emat * (Dzc * H) ;

H0 = 0.8 ; 
kappa = 1 ; 

fx = -kappa*(H - H0).*Hx ; 
fy = -kappa*(H - H0).*Hy ; 
fz = -kappa*(H - H0).*Hz ;


figure(1) ; clf ; 
Hplot = Eplot*H ;
Hplot = reshape(Hplot, size(xp));
surf(xp, yp, zp, Hplot, 'FaceAlpha', 0.5); hold on ;
shading interp
axis equal ;
title('numerical')



%% Inititalization of time-stepping

T = 4 ;

t = 0 ; 
dt = dx/5 ; 
numsteps = ceil(T/dt) ;
dt = T/numsteps ;

% y = [-1.209; 0; 0.341056] ; 
y = [-1.38826; 0.0436279; 0.0933781] ;

figure(1) ; plot3(y(1), y(2), y(3), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 7)

y_stored = zeros(3, numsteps) ;
count_proj =  0 ;


%% Time-stepping

for k = 1:numsteps

    
    %% Explicit Trapezoidal (RK2) time integration

    % Stage 1
    Etemp = interp3_matrix(x1d, y1d, z1d, y(1), y(2), y(3), 3, band) ;
    k1 = Etemp*[fx, fy, fz] ; k1 = k1' ;

    % Stage 2
    ytemp = y + dt*k1 ;
    Etemp = interp3_matrix(x1d, y1d, z1d, ytemp(1), ytemp(2), ytemp(3), 3, band) ; 
    k2 = Etemp*[fx, fy, fz] ; k2 = k2' ;

    % Update from TVD-RK2
    y = y + dt*(k1/2 + k2/2) ;
    
    %% Final Projection
    Etemp = interp3_matrix(x1d, y1d, z1d, y(1), y(2), y(3), 3, band) ;
    yproj = Etemp*[cpx_band, cpy_band, cpz_band] ; yproj = yproj' ;
    if norm(y - yproj) > dx/sqrt(2)
        y = yproj ; 
        count_proj = count_proj + 1 ; 
    end


    %% Save solution
    y_stored(:,k) = y ; 
   
end



figure(1)
plot3(y_stored(1,:), y_stored(2,:), y_stored(3,:), 'r-', 'LineWidth', 1.5)
axis equal 


% 
% % 
% figure(1) ; clf ; 
% Hplot = Eplot*H_exact ;
% Hplot = reshape(Hplot, size(xp));
% surf(xp, yp, zp, Hplot);
% axis equal ;
% title('exact')

figure(2) ; clf ; 
Eplot_test = interp3_matrix(x1d, y1d, z1d, vertices(:,1), vertices(:,2), vertices(:,3), 3, band) ;

trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3), Eplot_test*H, 'FaceAlpha', 1)
shading interp
axis equal