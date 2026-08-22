% Helix trajectory on torus

clear ; clc ; close all ;

jade = [0, 168, 107]/255 ;


%% Discretizaiton of embedding space

% Mesh size
dx = 0.1 ;

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
xx = xx(:) ; 
yy = yy(:) ; 
zz = zz(:) ;
cpxg = cpx(:) ;
cpyg = cpy(:) ;
cpzg = cpz(:) ;

% Computational band
dim = 3 ;
bw = 2*sqrt(dim) ;
band = find(abs(sdist) <= bw*dx) ;

% Store relevant closest point information in the band
xx_band = xx(band) ;
yy_band = yy(band) ; 
zz_band = zz(band) ;
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

y0 = y_exact(0) ;

%% Plot

figure(1) ; clf ;

% Plot surface
patch('Faces', face, 'Vertices', vertex, 'FaceVertexCData', ...
    0.5*ones(size(face)), 'EdgeColor', 'none', 'FaceColor', jade, 'FaceAlpha', 0.8) ;
hold on

% Plot exact solution on surface 
t_plot = linspace(0,2*pi,1000) ;
y_exact_plot = y_exact(t_plot) ;
plot3(y_exact_plot(1,:), y_exact_plot(2,:), y_exact_plot(3,:), '-k',...
    'Linewidth', 2) ;


axis equal
axis off

load('results/helix_torus/RK2_q=3_dx=0.04.mat') 
plot3(y_stored(1:2:end,1), y_stored(1:2:end,2), y_stored(1:2:end,3), 'ro', ...
      'MarkerSize', 5, 'MarkerFaceColor', 'r')

% Plot initial position
plot3(1.01*y0(1), 1.01*y0(2), 1.01*y0(3), 'bs', 'MarkerSize', 10, 'MarkerFaceColor', 'b')

view(-275.6672, 34.7296)

% exportgraphics(gcf, 'helix_torus_sol.pdf', 'Resolution', 600)