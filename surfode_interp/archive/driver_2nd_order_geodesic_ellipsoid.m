% Second-order geodesic equation on ellipsoid

clear ; clc ; close all ;


%% Discretization of embedding space

% Mesh size
dx = 0.025 ;

% Interpolation order
q = 1 ; 

% Make vectors of x, y, positions of the grid
x1d = (-1.2:dx:1.2)' ;
y1d = x1d ;
z1d = (-2.5:dx:2.5)' ;

% Embedding space
[xx, yy, zz] = meshgrid(x1d, y1d, z1d) ;

% Compute closest point 
R = 1 ;
r = 0.4 ;
[cpx, cpy, cpz, sdist] = cpEllipsoid(xx, yy, zz, [2 1], [], 'z') ;
[face, vertex] = isosurface(xx, yy, zz, sdist, 0) ;
xx = xx(:) ; 
yy = yy(:) ; 
zz = zz(:) ;
cpxg = cpx(:) ;
cpyg = cpy(:) ;
cpzg = cpz(:) ;
sdist = sdist(:) ;
% Computational band
dim = 3 ;
bw = sqrt(dim*((q + 1)/2)^2) ;
band = find(abs(sdist) <= bw*dx) ;

% Store relevant closest point information in the band
xx_band = xx(band) ;
yy_band = yy(band) ; 
zz_band = zz(band) ;
cpx_band = cpxg(band) ;
cpy_band = cpyg(band) ;
cpz_band = cpzg(band) ;
sdist_band = sdist(band) ;


%% Velocity field


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



[Dxx, Dyy, Dzz] = secondderiv_cen2_3d_matrices(x1d, y1d, z1d, band) ;
% 
[Dxy, Dxz, Dyz] = secondderiv_mixcen2_3d_matrices(x1d, y1d, z1d, band) ;
% 
% cpx_x = Dx*cpx_band ; 
% cpx_y = Dy*cpx_band ; 
% cpx_z = Dz*cpx_band ;
% 
% cpy_x = Dx*cpy_band ; 
% cpy_y = Dy*cpy_band ; 
% cpy_z = Dz*cpy_band ;

cpz_x = Dx*cpz_band ; 
cpz_y = Dy*cpz_band ; 
cpz_z = Dz*cpz_band ;

cpx_xx = Dxx*cpx_band ;
cpx_xy = Dxy*cpx_band ;
cpx_xz = Dxz*cpx_band ;
cpx_yy = Dyy*cpx_band ;
cpx_yz = Dyz*cpx_band ;
cpx_zz = Dzz*cpx_band ; 

cpy_xx = Dxx*cpy_band ;
cpy_xy = Dxy*cpy_band ;
cpy_xz = Dxz*cpy_band ;
cpy_yy = Dyy*cpy_band ;
cpy_yz = Dyz*cpy_band ;
cpy_zz = Dzz*cpy_band ; 

cpz_xx = Dxx*cpz_band ;
cpz_xy = Dxy*cpz_band ;
cpz_xz = Dxz*cpz_band ;
cpz_yy = Dyy*cpz_band ;
cpz_yz = Dyz*cpz_band ;
cpz_zz = Dzz*cpz_band ; 


%%

figure(1); clf;
patch('Faces', face, 'Vertices', vertex, 'FaceVertexCData', ...
    0.5*ones(size(face)), 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceAlpha', 0.5) ;
hold on;
axis equal ;
xlabel('x'); ylabel('y'); zlabel('z')

%% Time discretization


% Initialize time
t = 0 ;

% Simulation time
T = 20 ;

% Initial arrays for storing solution and error
y_stored = zeros(1,3) ;

dt = 0.001 ;

numtimesteps = ceil(T/dt) ;
x_stored = zeros(numtimesteps, 3) ;
u_stored = zeros(numtimesteps, 3) ;

% Initial condition
x = [0.7071067811865476 0.7071067811865475 1.2246467991473532e-16]' ;
u = [0.7020822129431468 0.7120871111210394 -0.007074586672431028]' - x ;

Etemp = interp3_matrix(x1d, y1d, z1d, x(1), x(2), x(3), q, band) ;
ntemp = [Etemp*nx; Etemp*ny; Etemp*nz] ; ntemp = ntemp/norm(ntemp) ;
uproj = u - dot(u, ntemp)*ntemp ;
u = uproj/norm(uproj) ;

% figure(1) ;
plot3(x(1), x(2), x(3), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 10)


% uproj = [Etemp*cpx_x Etemp*cpx_y Etemp*cpx_z ;
%          Etemp*cpy_x Etemp*cpy_y Etemp*cpy_z ;
%          Etemp*cpz_x Etemp*cpz_y Etemp*cpz_z ;]*u ;
% u = u/norm(u) ;


k1_x = zeros(3,1) ;
k1_u = zeros(3,1) ;
k2_x = zeros(3,1) ; 
k2_u = zeros(3,1) ;
k3_x = zeros(3,1) ; 
k3_u = zeros(3,1) ;
k4_x = zeros(3,1) ; 
k4_u = zeros(3,1) ;
count_proj = 0 ;



x0 = x ; 
du = zeros(3,1) ;
for k = 1:numtimesteps

    %% Forward Euler
    x = x0 + dt*u ;
    % Etemp = interp3_matrix(x1d, y1d, z1d, x0(1), x0(2), x0(3), q, band) ;
    % du(1) = u'*[Etemp*cpx_xx Etemp*cpx_xy Etemp*cpx_xz ;
    %             Etemp*cpx_xy Etemp*cpx_yy Etemp*cpx_yz ;
    %             Etemp*cpx_xz Etemp*cpx_yz Etemp*cpx_zz]*u ;
    % du(2) = u'*[Etemp*cpy_xx Etemp*cpy_xy Etemp*cpy_xz ;
    %             Etemp*cpy_xy Etemp*cpy_yy Etemp*cpy_yz ;
    %             Etemp*cpy_xz Etemp*cpy_yz Etemp*cpy_zz]*u ;
    % du(3) = u'*[Etemp*cpz_xx Etemp*cpz_xy Etemp*cpz_xz ;
    %             Etemp*cpz_xy Etemp*cpz_yy Etemp*cpz_yz ;
    %             Etemp*cpz_xz Etemp*cpz_yz Etemp*cpz_zz]*u ;
    % u = u + dt*du ;


    %% RK2
    % 
    % % Stage 1
    % k1_x = u ; 
    % Etemp = interp3_matrix(x1d, y1d, z1d, x(1), x(2), x(3), q, band) ;
    % k1_u(1) = u'*[Etemp*cpx_xx Etemp*cpx_xy Etemp*cpx_xz ;
    %               Etemp*cpx_xy Etemp*cpx_yy Etemp*cpx_yz ;
    %               Etemp*cpx_xz Etemp*cpx_yz Etemp*cpx_zz]*u ;
    % 
    % k1_u(2) = u'*[Etemp*cpy_xx Etemp*cpy_xy Etemp*cpy_xz ;
    %               Etemp*cpy_xy Etemp*cpy_yy Etemp*cpy_yz ;
    %               Etemp*cpy_xz Etemp*cpy_yz Etemp*cpy_zz]*u ;
    % 
    % k1_u(3) = u'*[Etemp*cpz_xx Etemp*cpz_xy Etemp*cpz_xz ;
    %               Etemp*cpz_xy Etemp*cpz_yy Etemp*cpz_yz ;
    %               Etemp*cpz_xz Etemp*cpz_yz Etemp*cpz_zz]*u ;
    % % k1_u = zeros(3,1) ;
    % 
    % 
    % % Stage 2
    % xtemp = x + dt*k1_x ; 
    % utemp = u + dt*k1_u ;    
    % Etemp = interp3_matrix(x1d, y1d, z1d, xtemp(1), xtemp(2), xtemp(3), q, band) ;
    % ntemp = [Etemp*nx; Etemp*ny; Etemp*nz] ; ntemp = ntemp/norm(ntemp) ;
    % utemp = utemp - dot(utemp, ntemp)*ntemp ;
    % k2_x = utemp ; 
    % k2_u(1) = utemp'*[Etemp*cpx_xx Etemp*cpx_xy Etemp*cpx_xz ;
    %                   Etemp*cpx_xy Etemp*cpx_yy Etemp*cpx_yz ;
    %                   Etemp*cpx_xz Etemp*cpx_yz Etemp*cpx_zz]*utemp ;
    % 
    % k2_u(2) = utemp'*[Etemp*cpy_xx Etemp*cpy_xy Etemp*cpy_xz ;
    %                   Etemp*cpy_xy Etemp*cpy_yy Etemp*cpy_yz ;
    %                   Etemp*cpy_xz Etemp*cpy_yz Etemp*cpy_zz]*utemp ;
    % 
    % k2_u(3) = utemp'*[Etemp*cpz_xx Etemp*cpz_xy Etemp*cpz_xz ;
    %                   Etemp*cpz_xy Etemp*cpz_yy Etemp*cpz_yz ;
    %                   Etemp*cpz_xz Etemp*cpz_yz Etemp*cpz_zz]*utemp ;
    % % k2_u = zeros(3,1) ;
    % 
    % % Update solution
    % x = x + dt*(k1_x/2 + k2_x/2) ; 
    % u = u + dt*(k1_u/2 + k2_u/2) ; 


    %% RK4

    % % Stage 1
    % k1_x = u ;
    % % k1_u(1) = u'*[Etemp*cpx_xx Etemp*cpx_xy Etemp*cpx_xz ;
    % %               Etemp*cpx_xy Etemp*cpx_yy Etemp*cpx_yz ;
    % %               Etemp*cpx_xz Etemp*cpx_yz Etemp*cpx_zz]*u ;
    % % 
    % % k1_u(2) = u'*[Etemp*cpy_xx Etemp*cpy_xy Etemp*cpy_xz ;
    % %               Etemp*cpy_xy Etemp*cpy_yy Etemp*cpy_yz ;
    % %               Etemp*cpy_xz Etemp*cpy_yz Etemp*cpy_zz]*u ;
    % % 
    % % k1_u(3) = u'*[Etemp*cpz_xx Etemp*cpz_xy Etemp*cpz_xz ;
    % %               Etemp*cpz_xy Etemp*cpz_yy Etemp*cpz_yz ;
    % %               Etemp*cpz_xz Etemp*cpz_yz Etemp*cpz_zz]*u ;
    % k1_u = zeros(3,1) ;
    % 
    % 
    % % Stage 2
    % xtemp = x + (dt/2)*k1_x ; 
    % utemp = u + (dt/2)*k1_u ;
    % Etemp = interp3_matrix(x1d, y1d, z1d, xtemp(1), xtemp(2), xtemp(3), q, band) ;
    % k2_x = utemp ; 
    % % k2_u(1) = utemp'*[Etemp*cpx_xx Etemp*cpx_xy Etemp*cpx_xz ;
    % %                   Etemp*cpx_xy Etemp*cpx_yy Etemp*cpx_yz ;
    % %                   Etemp*cpx_xz Etemp*cpx_yz Etemp*cpx_zz]*utemp ;
    % % 
    % % k2_u(2) = utemp'*[Etemp*cpy_xx Etemp*cpy_xy Etemp*cpy_xz ;
    % %                   Etemp*cpy_xy Etemp*cpy_yy Etemp*cpy_yz ;
    % %                   Etemp*cpy_xz Etemp*cpy_yz Etemp*cpy_zz]*utemp ;
    % % k2_u(3) = utemp'*[Etemp*cpz_xx Etemp*cpz_xy Etemp*cpz_xz ;
    % %                   Etemp*cpz_xy Etemp*cpz_yy Etemp*cpz_yz ;
    % %                   Etemp*cpz_xz Etemp*cpz_yz Etemp*cpz_zz]*utemp ;
    % k2_u = zeros(3,1) ;
    % 
    % 
    % % Stage 3
    % xtemp = x + (dt/2)*k2_x ; 
    % utemp = u + (dt/2)*k2_u ;
    % Etemp = interp3_matrix(x1d, y1d, z1d, xtemp(1), xtemp(2), xtemp(3), q, band) ;
    % k3_x = utemp ; 
    % % k3_u(1) = utemp'*[Etemp*cpx_xx Etemp*cpx_xy Etemp*cpx_xz ;
    % %                   Etemp*cpx_xy Etemp*cpx_yy Etemp*cpx_yz ;
    % %                   Etemp*cpx_xz Etemp*cpx_yz Etemp*cpx_zz]*utemp ;
    % % 
    % % k3_u(2) = utemp'*[Etemp*cpy_xx Etemp*cpy_xy Etemp*cpy_xz ;
    % %                   Etemp*cpy_xy Etemp*cpy_yy Etemp*cpy_yz ;
    % %                   Etemp*cpy_xz Etemp*cpy_yz Etemp*cpy_zz]*utemp ;
    % % k3_u(3) = utemp'*[Etemp*cpz_xx Etemp*cpz_xy Etemp*cpz_xz ;
    % %                   Etemp*cpz_xy Etemp*cpz_yy Etemp*cpz_yz ;
    % %                   Etemp*cpz_xz Etemp*cpz_yz Etemp*cpz_zz]*utemp ;
    % k3_u = zeros(3,1) ;
    % 
    % % Stage 4
    % xtemp = x + dt*k3_x ; 
    % utemp = u + dt*k3_u ;
    % Etemp = interp3_matrix(x1d, y1d, z1d, xtemp(1), xtemp(2), xtemp(3), q, band) ;
    % k4_x = utemp ; 
    % % k4_u(1) = utemp'*[Etemp*cpx_xx Etemp*cpx_xy Etemp*cpx_xz ;
    % %                   Etemp*cpx_xy Etemp*cpx_yy Etemp*cpx_yz ;
    % %                   Etemp*cpx_xz Etemp*cpx_yz Etemp*cpx_zz]*utemp ;
    % % 
    % % k4_u(2) = utemp'*[Etemp*cpy_xx Etemp*cpy_xy Etemp*cpy_xz ;
    % %                   Etemp*cpy_xy Etemp*cpy_yy Etemp*cpy_yz ;
    % %                   Etemp*cpy_xz Etemp*cpy_yz Etemp*cpy_zz]*utemp ;
    % % k4_u(3) = utemp'*[Etemp*cpz_xx Etemp*cpz_xy Etemp*cpz_xz ;
    % %                   Etemp*cpz_xy Etemp*cpz_yy Etemp*cpz_yz ;
    % %                   Etemp*cpz_xz Etemp*cpz_yz Etemp*cpz_zz]*utemp ;
    % k4_u = zeros(3,1) ;
    % 
    % 
    % % Update solution
    % x = x + dt*(k1_x/6 + k2_x/3 + k3_x/3 + k4_x/6) ; 
    % u = u + dt*(k1_u/6 + k2_u/3 + k3_u/3 + k4_u/6) ;


    %% Projection

    Etemp = interp3_matrix(x1d, y1d, z1d, x(1), x(2), x(3), q, band) ;
    xproj = [Etemp*cpx_band ; Etemp*cpy_band ; Etemp*cpz_band] ;
    ntemp = [Etemp*nx; Etemp*ny; Etemp*nz] ; ntemp = ntemp/norm(ntemp) ;
    u = u - dot(u, ntemp)*ntemp ; u = u/norm(u) ;
    x = xproj ;
    % if norm(x - xproj) > dx/sqrt(2) % || norm(u - uproj) > dx/sqrt(2)
    %     x = xproj ; 
    %     count_proj = count_proj + 1 ; 
    % end
    

    %% Save solution

    x_stored(k,:) = x' ;
    u_stored(k,:) = u' ;
    x0 = x ;

end


% disp(count_proj)
figure(1);
plot3(x_stored(:,1), x_stored(:,2), x_stored(:,3), 'r-', 'LineWidth', 2)
axis equal
xlabel('x'); ylabel('y'); zlabel('z')
save('colin/mytest.mat', 'x_stored')
