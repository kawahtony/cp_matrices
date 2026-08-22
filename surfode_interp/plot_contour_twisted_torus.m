clear ; clc ; close all ;

azure = [0, 128, 255]/255 ;
jade = [0, 168, 107]/255 ;
coral = [255, 127, 80]/255 ; 
indigo = [111, 0, 255]/255 ;

colorset = {azure, coral, indigo} ;
styleset = {'-', '--', ':'} ;


%% Load surface files

load('data/twisted_elliptic_torus_dx=0.025.mat') ;

xp = vertices(:,1) ;
yp = vertices(:,2) ;
zp = vertices(:,3) ;

Eplot = interp3_matrix(x1d, y1d, z1d, xp, yp, zp, q, band) ;
Hplot = Eplot*H ;


%% Plot mean curvature only

figure(1) ; clf ; hold on ;
trisurf(faces, xp, yp, zp, Hplot, 'FaceAlpha', 1, 'EdgeColor', 'none') ;
shading interp
colormap('cool')
view(273.4353, 38.3540)
axis equal off
camlight('headlight') ;
camlight('left') ;
lighting gouraud ;
material dull ;
colorbar('Location', 'southoutside', 'FontSize', 50)
exportgraphics(gcf, 'twisted_torus.pdf')



%% Plot contour trajectories with mean curvature map


% time_fraction = 1 ;  
% 
% figure(2) ; clf ; hold on ;
% 
% trisurf(faces, xp, yp, zp, Hplot, 'FaceAlpha', 0.25, 'EdgeColor', 'none') ;
% shading interp
% colormap('cool')
% 
% hLines = gobjects(1,3) ;
% for i = 1:3
%     load(['twisted_torus_contour_', num2str(i), '.mat'])
%     T = t_stored(end) ;
%     t_idx = find(t_stored >= time_fraction*T, 1) ;
%     hLines(i) = plot3(x_stored(1,1:t_idx), x_stored(2,1:t_idx), x_stored(3,1:t_idx), ...
%         styleset{i}, 'Color', colorset{i}, 'LineWidth', 4, 'DisplayName', num2str(i)) ;
%     plot3(x_stored(1,t_idx), x_stored(2,t_idx), x_stored(3,t_idx), 'o', ...
%         'MarkerFace', colorset{i}, 'MarkerSize', 7)
%     plot3(x_stored(1,1), x_stored(2,1), x_stored(3,1), 'o', ...
%         'MarkerEdge', 'k', 'LineWidth', 1.5, 'MarkerSize', 7)
% end
% 
% 
% % view(334.4864, 26.6017)
% view(273.4353, 38.3540)
% axis equal off
% fig_name = ['twisted_torus_fraction_', num2str(time_fraction), '.pdf'] ; 
% exportgraphics(gcf, fig_name)
