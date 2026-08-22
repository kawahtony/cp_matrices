% Test to find period of trajectories Euler's equations on sphere

clear ; clc ; close all


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

% Make vector
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


%% Looping over elevation angles

az = pi/3 ;
elset = linspace(pi/6, pi-pi/6, 5) ;

for k = 1:length(elset)

    % Initial location
    el = elset(k) ;
    y0 = [cos(az)*sin(el); sin(az)*sin(el); cos(el)] ;
    y = y0 ;
    
    % Intitalize time 
    t = 0 ;
    
    % Initialize CFL constant
    cfl_constant = max(f1) + max(f2) + max(f3) ;


    %% Intialize plot
    
    figure(1) ; clf ; 
    
    % Plot surface
    patch('Faces', face, 'Vertices', vertex, 'FaceVertexCData', ...
        0.5*ones(size(face)), 'EdgeColor', 'none', 'FaceColor', 'flat') ;
    hold on
    
    % Plot initial location
    plot3(y(1), y(2), y(3), 'ro', 'MarkerSize', 1)
    
    axis equal
    view(139,-5)
        
    
    %% Time-stepping
    
    k1 = zeros(3,1) ;
    k2 = zeros(3,1) ;
    
    res = 1 ;
    count = 0 ;
    
    while count <= 20 || res > 0.5*dx 
        

        %% Explicit Trapezoidal (RK2) time integrator
    
        dt = dx/cfl_constant ;
        
        % Stage 1 
        Etemp = interp3_matrix(x1d, y1d, z1d, y(1), y(2), y(3), 3, band) ;
        k1 = [Etemp*f1 ; Etemp*f2 ; Etemp*f3] ;
    
        % Stage 2
        ytemp = y + dt*k1 ;  
        Etemp = interp3_matrix(x1d, y1d, z1d, ytemp(1), ytemp(2), ytemp(3), 3, band) ;
        k2 = [Etemp*f1 ; Etemp*f2 ; Etemp*f3 ] ;
    
        % Update solution
        y = y + dt*(k1/2 + k2/2) ;
    
        % Projection
        Etemp = interp3_matrix(x1d, y1d, z1d, y(1), y(2), y(3), q, band) ;
        yproj = [Etemp*cpx_band ; Etemp*cpy_band ; Etemp*cpz_band] ;
        if norm(y - yproj) > dx/sqrt(2)
            y = yproj ; 
        end
    
        % Update cfl constant 
        cfl_constant = abs(Etemp*f1) + abs(Etemp*f2) + abs(Etemp*f3) ;
    
        % Update time
        t = t + dt ;

        % Update number of time steps
        count = count + 1 ;
    
        if count == 1
            y1 = y ;
        end

        % Update residual to the initial point
        res = min(norm(y - y0), norm(y - y1)) ;
        disp([t, res])
                
        % Plot numerical solution
        plot3(y(1), y(2), y(3), 'ro', 'MarkerSize', 3)
             
    
    end
    
    disp(t)

    
    save(['data/euler_sphere/traj_', num2str(k), '.mat'], 'y0', 't')
        
    pause

end


