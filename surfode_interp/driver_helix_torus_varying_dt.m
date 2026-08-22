% Helix trajectory on torus
% Version: Precompute embedding velocity in band

clear ; clc ; close all ;


%% Discretizaiton of embedding space

% Mesh size
dx = 0.01 ;

% Order of RK method (2, 3 or 4)
p = 4 ;

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
[face, vertex] = isosurface(xx, yy, zz, sdist, 0) ;
xx = xx(:) ; 
yy = yy(:) ; 
zz = zz(:) ;
cpxg = cpx(:) ;
cpyg = cpy(:) ;
cpzg = cpz(:) ;

% Computational band
dim = 3 ;
bw = 1 + sqrt(dim)*(q + 1)/2 ;
band = find(abs(sdist) <= bw*dx) ;

% Store relevant closest point information in the band
xx_band = xx(band) ;
yy_band = yy(band) ; 
zz_band = zz(band) ;
cpx_band = cpxg(band) ;
cpy_band = cpyg(band) ;
cpz_band = cpzg(band) ;


%% Velocity field

% Exact solution
omega1 = 6 ;
omega2 = 1 ;
x_exact = @(t) [(R + r*cos(omega1*t)).*cos(omega2*t); ...
                (R + r*cos(omega1*t)).*sin(omega2*t); ...
                r*sin(omega1*t)] ;


% Velocity field of helix trajectory
f1 = -omega1*zz_band.*cos(atan2(yy_band, xx_band)) - omega2*yy_band ;
f2 = -omega1*zz_band.*sin(atan2(yy_band, xx_band)) + omega2*xx_band ;
f3 = omega1*(xx_band.*cos(atan2(yy_band, xx_band)) + yy_band.*sin(atan2(yy_band,xx_band)) - R) ;

% cfl_constant = sqrt(2*(omega1*r + omega2*(R + r))^2 + (omega1*r)^2) ;
cfl_constant = (3*R + 4*r)*omega1 + 2*(R + r)*omega2 ; 


%% Time discretization


% Initialize time
t = 0 ;

% Simulation time
T = 2*pi ;

% Initial location
x = x_exact(0) ;

% Initial arrays for storing solution and error
x_stored = zeros(3,1) ;
x_exact_stored = zeros(3,1) ;
err_stored = zeros(1,1) ;
dt_stored = zeros(1,1) ;


%% Intialize plots


figure(1) ; clf ;

% Plot surface
patch('Faces', face, 'Vertices', vertex, 'FaceVertexCData', ...
    0.5*ones(size(face)), 'EdgeColor', 'none', 'FaceColor','flat') ;
hold on

% Plot exact solution on surface 
t_plot = linspace(0,2*pi,1000) ;
x_exact_plot = x_exact(t_plot) ;
plot3(x_exact_plot(1,:), x_exact_plot(2,:), x_exact_plot(3,:), '-k',...
    'Linewidth', 2) ;

% Plot initial position
plot3(x(1), x(2), x(3), 'bo', 'MarkerSize', 3)
axis equal
axis off


%% Time-stepping

k1 = zeros(3,1) ;
k2 = zeros(3,1) ;
k3 = zeros(3,1) ;
k4 = zeros(3,1) ;

count = 1 ; 
count_proj = 0 ; 

while t < T

    dt = dx/cfl_constant ;
        
    if dt > T-t
        dt = T-t ;
    end    

    switch p

        case 2 % Explicit Trapezoidal (RK2) time integrator            
            
            % Stage 1 
            Etemp = interp3_matrix(x1d, y1d, z1d, x(1), x(2), x(3), q, band) ;
            k1 = (Etemp*[f1, f2, f3])' ;
        
            % Stage 2
            ytemp = x + dt*k1 ;  
            Etemp = interp3_matrix(x1d, y1d, z1d, ytemp(1), ytemp(2), ytemp(3), q, band) ;
            k2 = (Etemp*[f1, f2, f3])' ;
        
            % Update solution
            x = x + dt*(k1/2 + k2/2) ;


        case 3 % Heun's 3rd order (RK3) time integrator 
            
            % Stage 1 
            Etemp = interp3_matrix(x1d, y1d, z1d, x(1), x(2), x(3), q, band) ;
            k1 = (Etemp*[f1, f2, f3])' ;
   
        
            % Stage 2 
            ytemp = x + (dt/3)*k1 ;
            Etemp = interp3_matrix(x1d, y1d, z1d, ytemp(1), ytemp(2), ytemp(3), q, band) ;
            k2 = (Etemp*[f1, f2, f3])' ;
        
            % Stage 3
            ytemp = x + (2*dt/3)*k2 ;
            Etemp = interp3_matrix(x1d, y1d, z1d, ytemp(1), ytemp(2), ytemp(3), q, band) ;
            k3 = (Etemp*[f1, f2, f3])' ;
        
            % Update of solution
            x = x + dt*(k1/4 + 3*k3/4) ;


        case 4 % Classical RK4 time integrator 
            
            % Stage 1 
            Etemp = interp3_matrix(x1d, y1d, z1d, x(1), x(2), x(3), q, band) ;
            k1 = (Etemp*[f1, f2, f3])' ;
        
            % Stage 2 
            ytemp = x + (dt/2)*k1 ;
            Etemp = interp3_matrix(x1d, y1d, z1d, ytemp(1), ytemp(2), ytemp(3), q, band) ;
            k2 = (Etemp*[f1, f2, f3])' ;
        
            % Stage 3
            ytemp = x + (dt/2)*k2 ; 
            Etemp = interp3_matrix(x1d, y1d, z1d, ytemp(1), ytemp(2), ytemp(3), q, band) ;
            k3 = (Etemp*[f1, f2, f3])' ;   
        
            % Stage 4
            ytemp = x + dt*k3 ; 
            Etemp = interp3_matrix(x1d, y1d, z1d, ytemp(1), ytemp(2), ytemp(3), q, band) ;
            k4 = (Etemp*[f1, f2, f3])' ;    
        
            % Update solution
            x = x + dt*(k1/6 + k2/3 + k3/3 + k4/6) ;

    end

    % Projection 
    Etemp = interp3_matrix(x1d, y1d, z1d, x(1), x(2), x(3), q, band) ;
    xproj = (Etemp*[cpx_band, cpy_band, cpz_band])' ;
    if norm(x - xproj) > dx/sqrt(2)
        x = xproj ; 
        count_proj = count_proj + 1 ; 
    end


    % Update cfl constant    
    cfl_constant = abs(Etemp*f1) + abs(Etemp*f2) + abs(Etemp*f3) ;


    % Plot numerical solution 
    plot3(x(1), x(2), x(3), 'ro', 'MarkerSize', 3, 'MarkerFaceColor', 'r')
    
    t = t + dt ;
    x_stored(:,count) = x ;
    x_exact_stored(:, count) = x_exact(t) ;
    dt_stored(count) = dt ;

    count = count + 1 ;
    
end

filename = ['results/helix_torus/RK', num2str(p), '_q=', num2str(q), '_dx=', num2str(dx), '.mat'] ;
save(filename, 'dx', 'x_stored', 'x_exact_stored', 'count_proj')