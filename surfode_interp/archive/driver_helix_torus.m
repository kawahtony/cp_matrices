% Helix trajectory on torus (interp for cp)

clear ; clc ; close all ;


%% Discretizaiton of embedding space

% Mesh size
dx = 0.0125 ;

% Make vectors of x, y, positions of the grid
x1d = (-2:dx:2)' ;
y1d = x1d ;
z1d = x1d ;
relpt = [x1d(1) y1d(1) z1d(1)] ;

% Embedding space
[xx, yy, zz] = meshgrid(x1d, y1d, z1d) ;

% Compute closest point funciton of torus
R = 1 ;
r = 0.4 ;
[cpx, cpy, cpz, sdist] = cpTorus(xx, yy, zz, R, r) ;
[face, vertex] = isosurface(xx, yy, zz, sdist, 0) ;
cpxg = cpx(:) ;
cpyg = cpy(:) ;
cpzg = cpz(:) ;

% Computational band
dim = 3 ;
bw = 2*sqrt(dim) ;
band = find(abs(sdist) <= bw*dx) ;

% Store relevant closest point information in the band
sdist_band = sdist(band) ; 
cpx_band = cpxg(band) ;
cpy_band = cpyg(band) ;
cpz_band = cpzg(band) ;


%% Time discretization

% Exact solution
omega1 = 6 ;
omega2 = 1 ;
y_exact = @(t) [(R + r*cos(omega1*t)).*cos(omega2*t); ...
                (R + r*cos(omega1*t)).*sin(omega2*t); ...
                r*sin(omega1*t)] ;


% Velocity field of helix trajectory
f = @(y) [-omega1*y(3).*cos(atan2(y(2),y(1))) - omega2*y(2);...
          -omega1*y(3).*sin(atan2(y(2),y(1))) + omega2*y(1);...
          omega1*(y(1).*cos(atan2(y(2),y(1))) + y(2).*sin(atan2(y(2),y(1))) - R)] ;
f_max = sqrt(2*(omega1*r + omega2*(R + r))^2 + (omega1*r)^2) ;

% Initialize time
t = 0 ;

% Simulation time
T = 2*pi ;

% Time step (CFL condition)
dt = dx/max(1,f_max) ;
numtimesteps = ceil(T/dt) ;
dt = T/numtimesteps ;

% Initial location of particles
y = y_exact(0) ;

% Initial arrays for storing solution and error
y_stored = zeros(numtimesteps,3) ;
y_exact_stored = zeros(numtimesteps,3) ;
err_stored = zeros(numtimesteps,1) ;


%% Intialize plots


figure(1) ; clf ;

% Plot surface
patch('Faces', face, 'Vertices', vertex, 'FaceVertexCData', ...
    0.5*ones(size(face)), 'EdgeColor', 'none', 'FaceColor','flat') ;
hold on

% Plot exact solution on surface 
t_plot = linspace(0,2*pi,1000) ;
y_exact_plot = y_exact(t_plot) ;
plot3(y_exact_plot(1,:), y_exact_plot(2,:), y_exact_plot(3,:), '-k',...
    'Linewidth', 2) ;

% Plot initial position
plot3(y(1), y(2), y(3), 'bo', 'MarkerSize', 3)
axis equal
axis off


%% Time-stepping

for kt = 1:numtimesteps

    %----------------------------------------------------------------------
    % Explicit Trapezoidal (RK2) time integrator
    
    % Stage 1 (assume y is on the surface)
    k1 = f(y) ;

    % Stage 2
    ytemp = y + dt*k1 ;  
    Etemp = interp3_matrix(x1d, y1d, z1d, ytemp(1), ytemp(2), ytemp(3), 3, band) ;
    ytemp(1) = Etemp*cpx_band ; 
    ytemp(2) = Etemp*cpy_band ;
    ytemp(3) = Etemp*cpz_band ;
    k2 = f(ytemp) ;

    % Update solution
    y = y + dt*(k1/2 + k2/2) ;

    % % Final projection 
    % Etemp = interp3_matrix(x1d, y1d, z1d, y(1), y(2), y(3), 3, band) ;
    % y(1) = Etemp*cpx_band ; 
    % y(2) = Etemp*cpy_band ;
    % y(3) = Etemp*cpz_band ;

    
    %----------------------------------------------------------------------
    % Heun's 3rd order (RK3) time integrator 
    
    % % Stage 1 (assume y is on the surface)
    % k1 = f(y) ;
    % 
    % % Stage 2 
    % ytemp = y + (dt/3)*k1 ;
    % Etemp = interp3_matrix(x1d, y1d, z1d, ytemp(1), ytemp(2), ytemp(3), 3, band) ;
    % ytemp(1) = Etemp*cpx_band ; 
    % ytemp(2) = Etemp*cpy_band ;
    % ytemp(3) = Etemp*cpz_band ;
    % k2 = f(ytemp) ;
    % 
    % % Stage 3
    % ytemp = y + (2*dt/3)*k2 ;
    % Etemp = interp3_matrix(x1d, y1d, z1d, ytemp(1), ytemp(2), ytemp(3), 3, band) ;
    % ytemp(1) = Etemp*cpx_band ; 
    % ytemp(2) = Etemp*cpy_band ;
    % ytemp(3) = Etemp*cpz_band ;
    % k3 = f(ytemp) ;
    % 
    % % Update of solution
    % y = y + dt*(k1/4 + 3*k3/4) ;
    % % 
    % % % Final projection 
    % % Etemp = interp3_matrix(x1d, y1d, z1d, y(1), y(2), y(3), 3, band) ;
    % % y(1) = Etemp*cpx_band ; 
    % % y(2) = Etemp*cpy_band ;
    % % y(3) = Etemp*cpz_band ;

    %----------------------------------------------------------------------
    % Classical RK4 time integrator 

    % % Stage 1 (assume y is on the surface)
    % k1 = f(y) ;
    % 
    % % Stage 2 
    % ytemp = y + (dt/2)*k1 ;
    % Etemp = interp3_matrix(x1d, y1d, z1d, ytemp(1), ytemp(2), ytemp(3), 3, band) ;
    % ytemp(1) = Etemp*cpx_band ; 
    % ytemp(2) = Etemp*cpy_band ;
    % ytemp(3) = Etemp*cpz_band ;
    % k2 = f(ytemp) ; 
    % 
    % % Stage 3
    % ytemp = y + (dt/2)*k2 ; 
    % Etemp = interp3_matrix(x1d, y1d, z1d, ytemp(1), ytemp(2), ytemp(3), 3, band) ;
    % ytemp(1) = Etemp*cpx_band ; 
    % ytemp(2) = Etemp*cpy_band ;
    % ytemp(3) = Etemp*cpz_band ;
    % k3 = f(ytemp) ; 
    % 
    % % Stage 4
    % ytemp = y + dt*k3 ; 
    % Etemp = interp3_matrix(x1d, y1d, z1d, ytemp(1), ytemp(2), ytemp(3), 3, band) ;
    % ytemp(1) = Etemp*cpx_band ; 
    % ytemp(2) = Etemp*cpy_band ;
    % ytemp(3) = Etemp*cpz_band ;
    % k4 = f(ytemp) ; 
    % 
    % % Update of solution
    % y = y + dt*(k1/6 + k2/3 + k3/3 + k4/6) ;
    % 
    % % % Final projection 
    % % Etemp = interp3_matrix(x1d, y1d, z1d, y(1), y(2), y(3), 3, band) ;
    % % y(1) = Etemp*cpx_band ; 
    % % y(2) = Etemp*cpy_band ;
    % % y(3) = Etemp*cpz_band ;


    %% Plot numerical solution 

    plot3(y(1), y(2), y(3), 'ro', 'MarkerSize', 3)
    t = kt*dt ;
    y_stored(kt, :) = y ;
    y_exact_stored(kt,:) = y_exact(t) ;
    
end

filename = ['results/helix_torus_RK2_dx=', num2str(dx), '.mat'] ;
save(filename, 'dt', 'dx', 'y_stored', 'y_exact_stored')
