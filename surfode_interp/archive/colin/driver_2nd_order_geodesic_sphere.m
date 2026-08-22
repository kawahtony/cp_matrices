% Second-order geodesic equation on ellipsoid

clear ; clc ; close all ;


%% Discretization of embedding space

% Mesh size
dx = 0.05 ;

% Interpolation order
q = 3 ; 

% Make vectors of x, y, positions of the grid
x1d = (-1.2:dx:1.2)' ;
y1d = x1d ;
z1d = x1d' ;

% Embedding space
[xx, yy, zz] = meshgrid(x1d, y1d, z1d) ;

% Compute closest point 
R = 1 ;
[cpx, cpy, cpz, sdist] = cpSphere(xx, yy, zz) ;
[face, vertex] = isosurface(xx, yy, zz, sdist, 0) ;
xx = xx(:) ; 
yy = yy(:) ; 
zz = zz(:) ;
cpxg = cpx(:) ;
cpyg = cpy(:) ;
cpzg = cpz(:) ;

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


%% Velocity field

nx = xx_band - cpx_band ; 
ny = yy_band - cpy_band ;
nz = zz_band - cpz_band ; 
nnorm = sqrt(nx.^2 + ny.^2 + nz.^2) ;
nx = nx./nnorm ; 
ny = ny./nnorm ; 
nz = nz./nnorm ;


[Dx, Dy, Dz] = firstderiv_cen2_3d_matrices(x1d, y1d, z1d, band) ;

[Dxx, Dyy, Dzz] = secondderiv_cen2_3d_matrices(x1d, y1d, z1d, band) ;

[Dxy, Dxz, Dyz] = secondderiv_mixcen2_3d_matrices(x1d, y1d, z1d, band) ;

cpx_x = Dx*cpx_band ; 
cpx_y = Dy*cpx_band ; 
cpx_z = Dz*cpx_band ;

cpy_x = Dx*cpy_band ; 
cpy_y = Dy*cpy_band ; 
cpy_z = Dz*cpy_band ;

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


%% Time discretization


% Initialize time
t = 0 ;

% Simulation time
T = 2*pi ;

% Initial arrays for storing solution and error
y_stored = zeros(1,3) ;

dt = 0.02 ;

numtimesteps = ceil(T/dt) ;
x_stored = zeros(numtimesteps, 3) ;
u_stored = zeros(numtimesteps, 3) ;

% Initial condition
x = [1/3; 2/3; 2/3] ;
u = [2; -1; 0]; 
u = u/norm(u) ;
figure(1)
plot3(x(1), x(2), x(3), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 20)

k1_x = zeros(3,1) ;
k1_u = zeros(3,1) ;
k2_x = zeros(3,1) ; 
k2_u = zeros(3,1) ;
k3_x = zeros(3,1) ; 
k3_u = zeros(3,1) ;
k4_x = zeros(3,1) ; 
k4_u = zeros(3,1) ;
count_proj = 0 ;

Etemp = interp3_matrix(x1d, y1d, z1d, x(1), x(2), x(3), q, band) ;

for k = 1:numtimesteps

    %% RK2

    % Stage 1
    k1_x = u ; 
    Etemp = interp3_matrix(x1d, y1d, z1d, x(1), x(2), x(3), q, band) ;
    k1_u(1) = u'*[Etemp*cpx_xx Etemp*cpx_xy Etemp*cpx_xz ;
                  Etemp*cpx_xy Etemp*cpx_yy Etemp*cpx_yz ;
                  Etemp*cpx_xz Etemp*cpx_yz Etemp*cpx_zz]*u ;

    k1_u(2) = u'*[Etemp*cpy_xx Etemp*cpy_xy Etemp*cpy_xz ;
                  Etemp*cpy_xy Etemp*cpy_yy Etemp*cpy_yz ;
                  Etemp*cpy_xz Etemp*cpy_yz Etemp*cpy_zz]*u ;

    k1_u(3) = u'*[Etemp*cpz_xx Etemp*cpz_xy Etemp*cpz_xz ;
                  Etemp*cpz_xy Etemp*cpz_yy Etemp*cpz_yz ;
                  Etemp*cpz_xz Etemp*cpz_yz Etemp*cpz_zz]*u ;
    % k1_u = zeros(3,1) ;


    % Stage 2
    xtemp = x + dt*k1_x ; 
    utemp = u + dt*k1_u ; 
    Etemp = interp3_matrix(x1d, y1d, z1d, xtemp(1), xtemp(2), xtemp(3), q, band) ;
    k2_x = utemp ; 
    k2_u(1) = utemp'*[Etemp*cpx_xx Etemp*cpx_xy Etemp*cpx_xz ;
                      Etemp*cpx_xy Etemp*cpx_yy Etemp*cpx_yz ;
                      Etemp*cpx_xz Etemp*cpx_yz Etemp*cpx_zz]*utemp ;

    k2_u(2) = utemp'*[Etemp*cpy_xx Etemp*cpy_xy Etemp*cpy_xz ;
                      Etemp*cpy_xy Etemp*cpy_yy Etemp*cpy_yz ;
                      Etemp*cpy_xz Etemp*cpy_yz Etemp*cpy_zz]*utemp ;

    k2_u(3) = utemp'*[Etemp*cpz_xx Etemp*cpz_xy Etemp*cpz_xz ;
                      Etemp*cpz_xy Etemp*cpz_yy Etemp*cpz_yz ;
                      Etemp*cpz_xz Etemp*cpz_yz Etemp*cpz_zz]*utemp ;
    % k2_u = zeros(3,1) ;

    % Update solution
    x = x + dt*(k1_x/2 + k2_x/2) ; 
    u = u + dt*(k1_u/2 + k2_u/2) ; 



    %% Projection

    Etemp = interp3_matrix(x1d, y1d, z1d, x(1), x(2), x(3), q, band) ;
    xproj = [Etemp*cpx_band ; Etemp*cpy_band ; Etemp*cpz_band] ;
    ntemp = x - xproj ; ntemp = ntemp/norm(ntemp) ;
    uproj = u - (ntemp'*u)*ntemp ;
    % uproj = [Etemp*cpx_x Etemp*cpx_y Etemp*cpx_z ;
    %          Etemp*cpy_x Etemp*cpy_y Etemp*cpy_z ;
    %          Etemp*cpz_x Etemp*cpz_y Etemp*cpz_z ;]*u ;
    x = xproj ;
    u = uproj ;
    % if norm(x - xproj) > dx/sqrt(2) % || norm(u - uproj) > dx/sqrt(2)
    %     x = xproj ; 
    %     % u = uproj ;
    %     count_proj = count_proj + 1 ; 
    % end
    u = u/norm(u) ;

    %% Update and save solution

    x_stored(k,:) = x' ;
    u_stored(k,:) = u' ;

    x0 = x ; 
    u0 = u ; 

end


disp(count_proj)
figure(1); clf;
patch('Faces', face, 'Vertices', vertex, 'FaceVertexCData', ...
      0.5*ones(size(face)), 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceAlpha', 0.5) ;
hold on;
plot3(x_stored(:,1), x_stored(:,2), x_stored(:,3), 'r-', 'LineWidth', 2)
axis equal
xlabel('x'); ylabel('y'); zlabel('z')
% save('colin/mytest.mat', 'x_stored')
