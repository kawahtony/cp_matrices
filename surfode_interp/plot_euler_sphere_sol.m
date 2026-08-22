% Euler's equations on sphere

clear ; clc ; close all

jade = [0, 168, 107]/255 ;


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
[cpx, cpy, cpz, sdist] = cpSphere(xx, yy, zz) ;
[face, vertex] = isosurface(xx, yy, zz, sdist, 0) ;


%% Plot

figure(1) ; clf ;

% Plot surface
patch('Faces', face, 'Vertices', vertex, 'FaceVertexCData', ...
    0.45*ones(size(face)), 'EdgeColor', 'none', 'FaceColor', jade, 'FaceAlpha', 0.8) ;
hold on


load('results/euler_sphere/traj_1.mat') 
plot3(x_stored(1,:), x_stored(2,:), x_stored(3,:), 'r-', 'LineWidth', 1)

load('results/euler_sphere/traj_2.mat') 
plot3(x_stored(1,:), x_stored(2,:), x_stored(3,:), 'r-', 'LineWidth', 1)

load('results/euler_sphere/traj_3.mat') 
plot3(x_stored(1,:), x_stored(2,:), x_stored(3,:), 'r-', 'LineWidth', 1)

load('results/euler_sphere/traj_4.mat') 
plot3(x_stored(1,:), x_stored(2,:), x_stored(3,:), 'r-', 'LineWidth', 1)

load('results/euler_sphere/traj_5.mat') 
plot3(x_stored(1,:), x_stored(2,:), x_stored(3,:), 'r-', 'LineWidth', 1)

axis equal
axis off

view(59.4187, 1.0918)

exportgraphics(gcf, 'euler_sphere.pdf')