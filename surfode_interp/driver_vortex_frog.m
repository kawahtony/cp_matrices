% %% Load obj 
% 
% 
% obj = readObj('data/frog.obj');
% vertices = obj.v;
% faces = obj.f.v;
% 
% figure(1); clf ;
% trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3), ...
%     'FaceColor', [0.4 0.7 0.4], 'EdgeColor', 'none');
% axis equal; camlight; lighting gouraud;
% title('Frog mesh');
% 
% 
% %% Construct a grid in the embedding space
% 
% % grid size
% dx = 0.025 ;
% 
% % make vectors of x, y, positions of the grid
% x1d = (-3:dx:3)' ;
% y1d = (-4:dx:3)' ;
% z1d = (-1:dx:4)' ;
% nx = length(x1d);
% ny = length(y1d);
% nz = length(z1d);
% 
% [xx, yy, zz] = meshgrid(x1d, y1d, z1d) ;
% 
% 
% %% Bandwidth
% 
% dim = 3 ; % Dimension of embedding space
% q = 3 ; % Degree of interpolating polynomial
% bw = 1 + sqrt(dim)*(q+1)/2 ; 
% 
% %%
% 
% disp('computing sdf')
% 
% [sdist, cpx, cpy, cpz, faceidx] = tri2sdf(xx, yy, zz, faces, vertices, bw*dx) ;
% 
% 
% disp('plotting')
% 
% figure(1) ; clf ;
% isosurface(xx, yy, zz, sdist, 0) ;
% axis equal ;
% 
% 
% 
% 
% %% Banding
% 
% cpxg = cpx(:) ;
% cpyg = cpy(:) ;
% cpzg = cpz(:) ;
% 
% band = find(abs(sdist) <= bw*dx) ;
% 
% % Store closest points in the band
% xx_band = xx(band) ;
% yy_band = yy(band) ;
% zz_band = zz(band) ;
% cpx_band = cpxg(band) ;
% cpy_band = cpyg(band) ;
% cpz_band = cpzg(band) ;
% 
% 
% save('data/frog_dx=0.025.mat')




%% Load sdf of frog surface directly

clear ; clc ; close all

jade = [0, 168, 107]/255 ;

load('data/frog_dx=0.025.mat')

[faces, vertices] = isosurface(xx, yy, zz, sdist, 0) ;

% figure(1) ; clf ;
% isosurface(xx, yy, zz, sdist, 0) ;
% axis equal ;


%% Velocity field

% f1 = -cpy_band ;
% f2 = cpx_band ; 
% f3 = cpz_band ;

f1 = 0 * cpy_band ;
f2 = cpz_band ; 
f3 = -cpx_band ;

% Projection on tangent space
for k = 1:length(band)
    vtemp = [f1(k); f2(k); f3(k)] ;
    ntemp = [nx(k); ny(k); nz(k)] ;
    vtemp = vtemp - dot(vtemp, ntemp)*ntemp ;
    f1(k) = vtemp(1) ;
    f2(k) = vtemp(2) ;
    f3(k) = vtemp(3) ;
end


%% Time discretization

% Initialize time 
t = 0 ;

% Initialize CFL condition
cfl_constant = max(abs(f1)) + max(abs(f2)) + max(abs(f3)) ;

% Simulation time
% T = 10 ;
T = 4 ;

% Initial location
x = [-1.125; -1.125; 2.41021] ;



%% Initialize plot

figure(1) ; clf ; 

% Plot surface
patch('Faces', faces, 'Vertices', vertices, 'FaceVertexCData', ...
    0.45*ones(size(faces)), 'EdgeColor', 'none', 'FaceColor', jade, 'FaceAlpha', 0.4) ;
hold on


% Plot initial location
plot3(x(1), x(2), x(3), 'ro', 'MarkerSize', 1)

axis equal
view(-48.3255, 22.2931)


%% Time stepping

k1 = zeros(3,1) ; 
k2 = zeros(3,1) ;
count_proj = 0 ; 
count = 1 ;
x_stored = x' ;


while t < T

    dt = dx / cfl_constant ;

    if T - t < dt 
        dt = T - t ;  
    end

    %% Explicit Trapezoidal (RK2) time integrator

    % Stage 1 
    Etemp = interp3_matrix(x1d, y1d, z1d, x(1), x(2), x(3), 3, band) ;
    k1 = [Etemp*f1 ; Etemp*f2 ; Etemp*f3] ;

    % Stage 2
    xtemp = x + dt*k1 ;  
    Etemp = interp3_matrix(x1d, y1d, z1d, xtemp(1), xtemp(2), xtemp(3), 3, band) ;
    k2 = [Etemp*f1 ; Etemp*f2 ; Etemp*f3 ] ;

    % Update solution
    x = x + dt*(k1/2 + k2/2) ;


    % Projection
    Etemp = interp3_matrix(x1d, y1d, z1d, x(1), x(2), x(3), q, band) ;
    xproj = [Etemp*cpx_band ; Etemp*cpy_band ; Etemp*cpz_band] ;


    % disp([t, y(1), y(2), y(3), yproj(1), yproj(2), yproj(3)])


    if norm(x - xproj) > dx/sqrt(2)
        x = xproj ; 
        count_proj = count_proj + 1 ;
    end

    % Update cfl constant 
    cfl_constant = abs(Etemp*f1) + abs(Etemp*f2) + abs(Etemp*f3) ;

    % Update number of time steps
    count = count + 1 ;

    % Update time
    t = t + dt ;

    % Save numerical solution
    x_stored(count, :) = x' ;


end

figure(1) ;
plot3(x_stored(:,1), x_stored(:,2), x_stored(:,3), 'r-', 'LineWidth', 1)


