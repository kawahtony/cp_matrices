% Euler's equations on sphere
    
clear ; clc ; close all


%% Discretization of embedding space

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
[cpx, cpy, cpz, sdist] = cpSphere(xx, yy, zz) ;
[face, vertex] = isosurface(xx, yy, zz, sdist, 0) ;

% Make vector
xx = xx(:) ; 
yy = yy(:) ; 
zz = zz(:) ;
cpxg = cpx(:) ;
cpyg = cpy(:) ;
cpzg = cpz(:) ;

% Computational band
dim = 3 ;
bw = 1 + sqrt(dim)*(q + 1)/2 ; ;
band = find(abs(sdist) <= bw*dx) ;

% Store closest points in the band
xx_band = xx(band) ;
yy_band = yy(band) ;
zz_band = zz(band) ;
cpx_band = cpxg(band) ;
cpy_band = cpyg(band) ;
cpz_band = cpzg(band) ;


%% Velocity field

% define velocity field
I1 = 1.6 ;
I2 = 1 ;
I3 = 2/3 ;

f1 = (1/I3 - 1/I2) * (yy_band.*zz_band) ;
f2 = (1/I1 - 1/I3) * (zz_band.*xx_band) ;
f3 = (1/I2 - 1/I1) * (xx_band.*yy_band) ;


%% Time discretization

traj_id = 5 ;

load(['data/euler_sphere/traj_', num2str(traj_id), '.mat']) ;
t_period = t ;

% Initialize time 
t = 0 ;

% Initialize CFL condition
cfl_constant = max(abs(f1)) + max(abs(f2)) + max(abs(f3)) ;

% Simulation time
T = 100 * t_period ;

% Initial location
% el = pi/4 ; az = pi/4 ;
% y = [cos(az)*sin(el); sin(az)*sin(el); cos(el)] ;
x = y0 ;

% Initialize count
count = 0 ;


%% Intialize plot

figure(1) ; clf ; 

% Plot surface
patch('Faces', face, 'Vertices', vertex, 'FaceVertexCData', ...
    0.5*ones(size(face)), 'EdgeColor', 'none', 'FaceColor', 'flat') ;
hold on

% Plot initial location
plot3(x(1), x(2), x(3), 'ro', 'MarkerSize', 1)

axis equal
view(139,-5)


%% Time-stepping

k1 = zeros(3,1) ;
k2 = zeros(3,1) ;
k3 = zeros(3,1) ;
k4 = zeros(3,1) ;

count_proj = 0 ; 

while t < T

    dt = dx/cfl_constant ;

    if T - t < dt
        dt = T - t ;
    end
        
    %% Explicit Trapezoidal (RK2) time integrator
    
    % Stage 1 
    Etemp = interp3_matrix(x1d, y1d, z1d, x(1), x(2), x(3), 3, band) ;
    k1 = [Etemp*f1 ; Etemp*f2 ; Etemp*f3] ;

    % Stage 2
    ytemp = x + dt*k1 ;  
    Etemp = interp3_matrix(x1d, y1d, z1d, ytemp(1), ytemp(2), ytemp(3), 3, band) ;
    k2 = [Etemp*f1 ; Etemp*f2 ; Etemp*f3 ] ;

    % Update solution
    x = x + dt*(k1/2 + k2/2) ;

    % Projection
    Etemp = interp3_matrix(x1d, y1d, z1d, x(1), x(2), x(3), q, band) ;
    yproj = [Etemp*cpx_band ; Etemp*cpy_band ; Etemp*cpz_band] ;
    if norm(x - yproj) > dx/sqrt(2)
        x = yproj ; 
        count_proj = count_proj + 1 ;
    end

    % Update cfl constant 
    cfl_constant = abs(Etemp*f1) + abs(Etemp*f2) + abs(Etemp*f3) ;
    
    % Update time
    t = t + dt ;

    % Update number of time steps
    count = count + 1 ;

    % Save numerical solution
    x_stored(:, count) = x ;
   

end


plot3(x_stored(1,:), x_stored(2,:), x_stored(3,:), 'r-', 'LineWidth', 1)


save(['results/euler_sphere/traj_', num2str(traj_id), '.mat'], 'dx', 'x_stored', 'count_proj')