% Gradient flow to trace geodesic path from geodesic field on bunny

clear ; clc ; close all ;


%% Load geodesic field

load('data/bunny_dx=0.025.mat')


%% Velocity field

phix = Eqmat*(Dxc * phi) ; 
phiy = Eqmat*(Dyc * phi) ;
phiz = Eqmat*(Dzc * phi) ;

fx = -phix ; 
fy = -phiy ;
fz = -phiz ; 


%% Inititalization of time-stepping

T = 2 ;

t = 0 ; 
dt = dx/2 ; 
numsteps = ceil(T/dt) ;
dt = T/numsteps ;

% y = [0.0231704; -0.240589; 0.857734] ; % ear
% y = [0.713078; 0.356921; -0.388132] ; % butt
% y = [-0.283861; 0.51601; -0.780162] ; % foot
y = [-0.837554; 0.211496; -0.321818] ; % side

y_stored = zeros(3, numsteps) ;
count_proj =  0 ;

for k = 1:numsteps

    
    %% Explicit Trapezoidal (RK2) time integration

    % Stage 1
    Etemp = interp3_matrix(x1d, y1d, z1d, y(1), y(2), y(3), 3, band) ;
    k1 = Etemp*[fx, fy, fz] ; k1 = k1' ;

    % Stage 2
    ytemp = y + dt*k1 ;
    Etemp = interp3_matrix(x1d, y1d, z1d, ytemp(1), ytemp(2), ytemp(3), 3, band) ; 
    k2 = Etemp*[fx, fy, fz] ; k2 = k2' ;

    % Update from TVD-RK2
    y = y + dt*(k1/2 + k2/2) ;
    
    %% Final Projection
    Etemp = interp3_matrix(x1d, y1d, z1d, y(1), y(2), y(3), 3, band) ;
    yproj = Etemp*[cpx_band, cpy_band, cpz_band] ; yproj = yproj' ;
    if norm(y - yproj) > dx/sqrt(2)
        y = yproj ; 
        count_proj = count_proj + 1 ; 
    end


    %% Save solution
    y_stored(:,k) = y ; 
   
end

figure(1) ; clf ; 
trisurf(Faces,xp,yp,zp, Eplot*phi, 'FaceAlpha', 0.5); hold on ;
shading interp
plot3(y_stored(1,:), y_stored(2,:), y_stored(3,:), 'r-', 'LineWidth', 1.5)
% quiver3(xp, yp, zp, X_x_plot, X_y_plot, X_z_plot, 'Color', 'r')
axis equal 

count = 1 ;
filename = ['results/geodesic_bunny_', num2str(count), '.mat'] ;
while exist(filename, 'file')
    count = count + 1 ;
    filename = ['results/geodesic_bunny_', num2str(count), '.mat'] ;
end

save(filename, 'y_stored', 'count_proj') 