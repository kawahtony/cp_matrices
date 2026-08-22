clear ; clc ; close all ;

jade = [0, 168, 107]/255 ;


%% Load surface files

load('data/wavy_torus_dx=0.025.mat') ;


%% Compute mean curvature (use unit outward normal)

Lmat = laplacian_3d_matrix(x1d, y1d, z1d, 2, band);
Emat = interp3_matrix(x1d, y1d, z1d, cpx_band, cpy_band, cpz_band, q, band);
[Dxc, Dyc, Dzc] = firstderiv_cen2_3d_matrices(x1d, y1d, z1d, band) ;

Lcpx = Emat * (Lmat * cpx_band) ; 
Lcpy = Emat * (Lmat * cpy_band) ; 
Lcpz = Emat * (Lmat * cpz_band) ; 

H = zeros(size(band)) ;

for k = 1:length(band)
    H(k) = -dot([Lcpx(k); Lcpy(k); Lcpz(k)], [nx(k); ny(k); nz(k)]) ;
end

H = Emat * H ; 


%% Initialize plot

xp = vertices(:,1) ;
yp = vertices(:,2) ;
zp = vertices(:,3) ;

Eplot = interp3_matrix(x1d, y1d, z1d, xp, yp, zp, q, band) ;
Hplot = Eplot*H ;


%% Plot surface only

scale = 0.985 ;

% figure(1) ; clf ; 
% patch('Faces', faces, 'Vertices', scale*vertices, 'FaceVertexCData', ...
%     0.45*ones(size(faces)), 'EdgeColor', 'none', 'FaceColor', jade, 'FaceAlpha', 1) ;
% axis equal off  
% camlight('headlight') ;
% camlight('left') ;
% lighting gouraud ;
% material dull ;
% view(0.2000, 58.3811) 
% exportgraphics(gcf, 'wavy_torus.pdf')




%% Critical points: ascent (maxima) and descent (minima)
% Each point is drawn with a white "halo" behind it so it pops against
% the curvature colormap regardless of the local surface color.


figure(2) ; clf ; hold on ;
trisurf(faces, scale*xp, scale*yp, scale*zp, Hplot, 'FaceAlpha', 1, 'EdgeColor', 'none') ;
shading interp
axis equal off  
cb = colorbar ;
set(cb, 'FontSize', 20)
cb.Label.Interpreter = 'Latex' ;
colormap cool
camlight('headlight') ;
camlight('left') ;
lighting gouraud ;
material dull ;
view(0.2000, 58.3811) 


load('critical_temp_ascent.mat')
xa = critical_stored(1,:) ; ya = critical_stored(2,:) ; za = critical_stored(3,:) ;
plot3(xa, ya, za, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'w', ...
    'MarkerSize', 18) ; % halo
h1 = plot3(xa, ya, za, 'o', 'MarkerFaceColor', [0.85, 0, 0], 'MarkerEdgeColor', 'k', ...
    'MarkerSize', 14, 'LineWidth', 2) ;

load('critical_temp_descent.mat')
xd = critical_stored(1,:) ; yd = critical_stored(2,:) ; zd = critical_stored(3,:) ;
plot3(xd, yd, zd, 's', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'w', ...
    'MarkerSize', 18) ; % halo
h2 = plot3(xd, yd, zd, 's', 'MarkerFaceColor', [0.1, 0.3, 0.9], 'MarkerEdgeColor', 'k', ...
    'MarkerSize', 14, 'LineWidth', 2) ;



legend([h1, h2], {'Ascent', 'Descent'}, 'Location', 'best', 'FontSize', 20, 'Box', 'off') ;

exportgraphics(gcf, 'wavy_torus_2.pdf')





