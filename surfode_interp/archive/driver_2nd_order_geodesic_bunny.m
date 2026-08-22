% Gradient flow to trace geodesic path from geodesic field on bunny

clear ; clc ; close all ;


%% Construct a grid in the embedding space

% grid size
dx = 0.05 ;

% make vectors of x, y, positions of the grid
x1d = (-2:dx:2)' ;
y1d = x1d ;
z1d = x1d ;

[xx, yy, zz] = meshgrid(x1d, y1d, z1d) ;


%% Bandwidth

dim = 3 ; % Dimension of embedding space
q = 3 ; % Degree of interpolating polynomial
fd_stenrad = 1 ; % Finite difference stencil radius (meaningless here)
bw = 1.00001*sqrt((dim-1)*((q+1)/2)^2 + ((fd_stenrad+(q+1)/2)^2)) ;
bw = bw * dx;


%% Find closest points on the surface

PlyFile = 'bunny.ply';

ptCloud = pcread(PlyFile) ;
Vertices = double(ptCloud.Location) ;
temp = Vertices(:,2) ;
Vertices(:,2) = Vertices(:,3) ; 
Vertices(:,3) = temp ; 
clear temp ;

Faces = MyRobustCrust(Vertices) ;
Faces = double(Faces) ;

% Compute closest point
[sdist, cpx, cpy, cpz, faceidx] = tri2sdf(xx, yy, zz, Faces, Vertices, bw) ;


figure(1) ; clf ; hold on ;
p = patch(isosurface(xx, yy, zz, sdist, 0)) ;
set(p, 'FaceColor', [0.85 0.82 0.78], 'EdgeColor', 'none') ;
isonormals(xx, yy, zz, sdist, p) ;
axis equal ;
view(-148.6109, 15.2982)
camlight('headlight') ; lighting gouraud ; 
material dull ;


% Reshape to vector
xx = xx(:) ; 
yy = yy(:) ;
zz = zz(:) ; 
cpx = cpx(:) ;
cpy = cpy(:) ;
cpz = cpz(:) ;
sdist = sdist(:) ;

% Computational band
band = find(abs(sdist) <= bw) ;

% Closest point information in the band
xx_band = xx(band) ;
yy_band = yy(band) ; 
zz_band = zz(band) ;
cpx_band = cpx(band) ;
cpy_band = cpy(band) ;
cpz_band = cpz(band) ;
sdist_band = sdist(band) ;


%% Compute unit outward normal

[Dx, Dy, Dz] = firstderiv_cen2_3d_matrices(x1d, y1d, z1d, band) ;

nx = (xx_band - cpx_band).*sign(sdist_band) ; 
ny = (yy_band - cpy_band).*sign(sdist_band)  ;
nz = (zz_band - cpz_band).*sign(sdist_band)  ; 
nnorm = sqrt(nx.^2 + ny.^2 + nz.^2) ;

phi_x = Dx*sdist_band ;
phi_y = Dy*sdist_band ;
phi_z = Dz*sdist_band ;
for k = 1:length(nnorm)
    if nnorm(k) < 1e-10
        nx(k) = phi_x(k) ; 
        ny(k) = phi_y(k) ;
        nz(k) = phi_z(k) ;
        nnorm(k) = sqrt(nx(k)^2 + ny(k)^2 + nz(k)^2) ;
    end
end

nx = nx./nnorm ; 
ny = ny./nnorm ; 
nz = nz./nnorm ;

% figure(1) ;
% quiver3(cpx_band, cpy_band, cpz_band, nx, ny, nz, 5, 'LineWidth', 0.1)


%% Time discretization 

% Simulation time 
T = 1 ;

% Time step 
dt = 0.00005 ;

numtimesteps = ceil(T/dt) ;
dt = T/numtimesteps ;
x_stored = zeros(numtimesteps, 3) ;

% Load a geodesic path (gradient flow from CPHM solution)
load('results/geodesic_bunny_3.mat', 'y_stored')
x = y_stored(:,1) ;
figure(1) ;
plot3(x(1), x(2), x(3), 'ro', 'MarkerFaceColor', 'r')
u = y_stored(:,2) - y_stored(:,1) ;
Etemp = interp3_matrix(x1d, y1d, z1d, x(1), x(2), x(3), q, band) ;
ntemp = [Etemp*nx; Etemp*ny; Etemp*nz] ; ntemp = ntemp/norm(ntemp) ;
u = u - dot(u, ntemp)*ntemp ; u = u/norm(u) ;


count_proj = 0 ;
t = 0 ;

tic ; 
for k = 1:numtimesteps

    % % Forward Euler timestepping
    x = x + dt*u ;

    % % Explicit Midpoint (RK2)
    % xtemp = x + (dt/2)*u ; 
    % Etemp = interp3_matrix(x1d, y1d, z1d, xtemp(1), xtemp(2), xtemp(3), q, band) ;
    % ntemp = [Etemp*nx; Etemp*ny; Etemp*nz] ; ntemp = ntemp/norm(ntemp) ;
    % u = u - dot(u, ntemp)*ntemp ; u = u/norm(u) ;
    % x = x + dt*u ;

    % Projection
    Etemp = interp3_matrix(x1d, y1d, z1d, x(1), x(2), x(3), q, band) ;
    xproj = [Etemp*cpx_band ; Etemp*cpy_band ; Etemp*cpz_band] ;
    ntemp = [Etemp*nx; Etemp*ny; Etemp*nz] ; ntemp = ntemp/norm(ntemp) ;
    u = u - dot(u, ntemp)*ntemp ; u = u/norm(u) ;
    if norm(x - xproj) > dx/sqrt(2) 
        x = xproj ; 
        count_proj = count_proj + 1 ; 
    end
    t = t + dt ;

    % Save solution
    x_stored(k,:) = x' ;

end
toc ; 

figure(1);
plot3(x_stored(:,1), x_stored(:,2), x_stored(:,3), 'r--', 'LineWidth', 2)
plot3(y_stored(1,:), y_stored(2,:), y_stored(3,:), 'k:', 'LineWidth', 1.5)


% Comparison between CPHM gradient flow solution (sol1) and second-order
% IVP solution (sol2) is not well-defined. Sol2 uses a normal computed from
% triangular face. Sol1 corresponds an application in geometry processing -
% retrieve geodesic path from geodesic distance field