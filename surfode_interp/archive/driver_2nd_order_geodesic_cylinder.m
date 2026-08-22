% Second-order geodesic equation on ellipsoid

clear ; clc ; close all ;
format long

%% Discretization of embedding space

% Mesh size
dx = 0.05 ;

% Interpolation order
q = 3 ; 

% Make vectors of x, y, positions of the grid
x1d = (-1.2:dx:1.2)' ;
y1d = x1d ;
z1d = (-2.3:dx:2.3)' ;

% Embedding space
[xx, yy, zz] = meshgrid(x1d, y1d, z1d) ;

% Compute closest point 
R = 0.7 ;
r = 0.4 ;
[cpx, cpy, cpz, dist] = cpCylinder(xx, yy, zz, [-2 2], R) ;
xx = xx(:) ; 
yy = yy(:) ; 
zz = zz(:) ;
cpxg = cpx(:) ;
cpyg = cpy(:) ;
cpzg = cpz(:) ;
dist = dist(:) ;

% Computational band
dim = 3 ;
bw = sqrt(dim*((q + 1)/2)^2) ;
band = find(dist <= bw*dx) ;

% Store relevant closest point information in the band
xx_band = xx(band) ;
yy_band = yy(band) ; 
zz_band = zz(band) ;
cpx_band = cpxg(band) ;
cpy_band = cpyg(band) ;
cpz_band = cpzg(band) ;
dist_band = dist(band) ;
sdist_band = dist_band .* sign(sqrt(xx_band.^2 + yy_band.^2) - R) ;


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


%%

[xp, yp, zp] = paramCylinder(200, [-2 2], 0.7) ;

figure(1); clf;
surf(xp, yp, zp, 'FaceAlpha', 0.5, 'EdgeColor', 'None')
hold on;
axis equal ;
xlabel('x'); ylabel('y'); zlabel('z')
colormap gray


%% Time discretization


c = 0.2 ; 
omega = sqrt((1 - c^2)/R^2) ; 

% Simulation time
T = 15 ;

t = linspace(0, T, 300)' ;
thetat = pi/3 + omega*t/(R^2*omega^2 + c^2) ;
zt = -1.5 + c*t/(R^2*omega^2 + c^2) ; 
xt = R*cos(thetat) ; 
yt = R*sin(thetat) ;

plot3(xt, yt, zt, 'k-', 'LineWidth', 2)


% Time step
dt = 0.01 ;

numtimesteps = ceil(T/dt) ;
dt = T/numtimesteps ;
x_stored = zeros(numtimesteps, 3) ;

% Initial condition
x = [xt(1); yt(1); zt(1)] ;
u = [-R*omega*sin(thetat(1)); R*omega*cos(thetat(1)); c] ; u = u/norm(u) ;
plot3(x(1), x(2), x(3), 'ro', 'MarkerFaceColor', 'r')

count_proj = 0 ;
t = 0 ;


for k = 1:numtimesteps

    % Forward Euler timestepping
    x = x + dt*u ;
    
    % % Explicit Midpoint (RK2)
    % xtemp = x + (dt/2)*u ; 
    % Etemp = interp3_matrix(x1d, y1d, z1d, xtemp(1), xtemp(2), xtemp(3), q, band) ;
    % ntemp = [Etemp*nx; Etemp*ny; Etemp*nz] ; ntemp = ntemp/norm(ntemp) ;
    % u = u - dot(u, ntemp)*ntemp ; u = u/norm(u) ;
    % x = x + dt*u ;

    % % Explicit Trapezoidal (RK2) 
    % k1 = u ;
    % xtemp = x + dt*k1 ;
    % Etemp = interp3_matrix(x1d, y1d, z1d, xtemp(1), xtemp(2), xtemp(3), q, band) ;
    % ntemp = [Etemp*nx; Etemp*ny; Etemp*nz] ; ntemp = ntemp/norm(ntemp) ;
    % k2 = u - dot(u, ntemp)*ntemp ; k2 = k2/norm(k2) ;
    % x = x + dt*(k1/2 + k2/2) ;
    
   
    % Projection
    Etemp = interp3_matrix(x1d, y1d, z1d, x(1), x(2), x(3), q, band) ;
    xproj = [Etemp*cpx_band ; Etemp*cpy_band ; Etemp*cpz_band] ;
    x = xproj ;
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

figure(1);
plot3(x_stored(:,1), x_stored(:,2), x_stored(:,3), 'r--', 'LineWidth', 2)
plot3(x(1), x(2), x(3), 'ro', 'MarkerFaceColor', 'r')
plot3(xt(end), yt(end), zt(end), 'ks', 'MarkerFaceColor', 'k')
disp(norm(x - [xt(end); yt(end); zt(end)]))



% k1_x = u ;
% k1_u = zeros(3,1) ;
% 
% xtemp = x + dt*k1_x ;
% utemp = u + dt*k1_u ; 
% Etemp = interp3_matrix(x1d, y1d, z1d, xtemp(1), xtemp(2), xtemp(3), q, band) ;
% ntemp = [Etemp*nx; Etemp*ny; Etemp*nz] ; ntemp = ntemp/norm(ntemp) ;
% k2_x = utemp - dot(utemp, ntemp)*ntemp ;
% 
% x = x + dt*(k1_x/2 + k2_x/2) ;