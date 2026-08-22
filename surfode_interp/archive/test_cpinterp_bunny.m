% Gradient flow to trace geodesic path from geodesic field on bunny

clear ; clc ; close all ;


%% Load geodesic field

load('data/bunny_dx=0.025.mat')


%% Velocity field

numtimesteps = 3000 ;
dt = T/numtimesteps ;

% y = [0.0231704; -0.240589; 0.857734] ; % ear
% y = [0.713078; 0.356921; -0.388132] ; % butt
% y = [-0.283861; 0.51601; -0.780162] ; % foot
y = [-0.837554; 0.211496; -0.321818] ; % side

y_stored = zeros(3, numtimesteps) ;


for k = 1:numtimesteps

    
    %% Projection
    Etemp = interp3_matrix(x1d, y1d, z1d, y(1), y(2), y(3), 3, band) ;
    y = Etemp*[cpx_band, cpy_band, cpz_band] ; y = y' ;


    %% Save solution
    y_stored(:,k) = y ; 
   
end



figure(1) ; clf ; 
trisurf(Faces,xp,yp,zp, Eplot*phi, 'FaceAlpha', 0.5); hold on ;
shading interp
plot3(y_stored(1,:), y_stored(2,:), y_stored(3,:), 'r-', 'LineWidth', 1.5)
axis equal 
