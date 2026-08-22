
%% n-body on trefoil knot tube

clc ; clear ; close all ; 

rng(42)

jade = [0, 168, 107]/255 ; 


%% Load closest point of trefoil knot tubulbar surface

load('data/trefoil_knot_tube_dx=0.025.mat')


%% Set-up for n-body simulation

% define velocity field
f = @velocity_field_nbody ;

% number of particles
num = 20 ;

% time-stepping initialization
t = 0 ;
T = 100 ;
x0_ind = randperm(length(band), num) ;
x = [cpx_band(x0_ind)'; cpy_band(x0_ind)'; cpz_band(x0_ind)'] ;
v = randn(3,num) ;


%%


% % setup for solution storage
% file_count = 1 ;
% folderName = strcat(shape, '_n=', num2str(num), ...
%     '_T=', num2str(T), '_example', num2str(file_count)) ;
% 
% while exist(folderName, 'file') == 7
%     file_count = file_count + 1 ;
%     folderName =  strcat(shape, '_n=', ...
%         num2str(num), '_T=', num2str(T), '_example', num2str(file_count)) ;    
% end
% 
% mkdir(folderName)
% savdir = strcat(pwd, '/', folderName) ;
dt_record = 0.01 ;
t_record = dt_record ;
count_record = 0 ;


%% Initialize plot

xp = vertices(:,1) ;
yp = vertices(:,2) ;
zp = vertices(:,3) ;

Eplot = interp3_matrix(x1d, y1d, z1d, xp, yp, zp, q, band) ;

figure(1) ; clf ; hold on ;
trisurf(faces, xp, yp, zp, Eplot*H, 'FaceAlpha', 1) ;
shading interp
axis equal tight off 
view(3)


% figure(1) ; clf ; hold on ;
% p = patch('Faces', faces, 'Vertices', vertices, 'FaceAlpha', 0.5, 'FaceColor', jade, ...
%           'EdgeColor', 'none') ;
% quiver3(cpx_band, cpy_band, cpz_band, nx, ny ,nz)
% axis equal tight off 

% % distinct color per particle, shared by its marker and its trail
% colors = lines(num) ;
% 
% % one animatedline per particle for the trail history: addpoints()
% % appends in place, so the trail never needs to be deleted/redrawn from
% % scratch as it grows
% trailMaxPts = 50 ;   % cap trail length so draw cost stays bounded over long runs
% trails = gobjects(num,1) ;
% for k = 1:num
%     trails(k) = animatedline('Color', colors(k,:), 'LineWidth', 1, ...
%         'MaximumNumPoints', trailMaxPts) ;
%     addpoints(trails(k), x(1,k), x(2,k), x(3,k)) ;
% end
% 
% % current particle positions; updated via set(...,'XData',...) each frame
% % instead of delete()+replot
% m1 = scatter3(x(1,:), x(2,:), x(3,:), 36, colors, 'filled', 'MarkerEdgeColor', 'k') ;
% 
% titleHandle = title('t = 0', 'FontSize', 16) ;
% 
% % view(-181.5807, 17.7373)
% 
% k1_x = zeros(size(x)) ;
% k1_v = zeros(size(x)) ;
% 
% alpha = 1 ; 
% beta = 0.5 ;
% 
% while t < T
% 
%     fmax = 1 ;
%     for k = 1:num
%         fmax = max(fmax, norm(f(x,k))) ;
%     end  
% 
%     dt = dx/fmax ; disp(dt)
% 
%     if T - t < dt
%         dt = T - t ;
%     end
% 
%     % stage 1
% 
%     for k = 1:num      
%         Etemp = interp3_matrix(x1d, y1d, z1d, x(1,k), x(2,k), x(3,k), 3, band) ;
%         ntemp = ( Etemp*[nx, ny, nz] )' ; 
%         k1_x(:,k)= v(:,k) ;
%         k1_v(:,k) = f(x,k) ;
%         k1_v(:,k) = k1_v(:,k) - dot(k1_v(:,k), ntemp)*ntemp ;
% 
%     end
%     % % stage 2
%     % xtemp = x + dt*k1 ;
%     % k2 = zeros(size(x)) ;
%     % for k = 1:num              
%     %     Etemp = interp3_matrix(x1d, y1d, z1d, xtemp(1,k), xtemp(2,k), xtemp(3,k), 3, band) ;
%     %     ntemp = Etemp*[nx, ny, nz] ; ntemp = ntemp' ; 
%     %     k2(:,k) = f(xtemp,k) ;
%     %     k2(:,k) = k2(:,k) - dot(k2(:,k), ntemp)*ntemp ;             
%     % end
%     % 
%     % % Update from RK2
%     % x = x + dt*(k1/2 + k2/2) ;
% 
%     % Update from FE
%     x = x + dt*k1_x ;
%     v = v + dt*k1_v ;
% 
%     % Projection
%     for k = 1:num
%         Etemp =  interp3_matrix(x1d, y1d, z1d, x(1,k), x(2,k), x(3,k), 3, band) ;
%         ntemp = ( Etemp*[nx, ny, nz] )' ;
%         xproj = ( Etemp*[cpx_band, cpy_band, cpz_band] )' ;
%         x(:,k) = xproj ;
%         v(:,k) = v(:,k) - dot(v(:,k), ntemp)*ntemp ;
%         % if norm(x(:,k) - xproj) > (0.1*bw)*dx
%         %     x(:,k) = xproj ;
%         % end
%     end
%     t = t + dt ;
% 
%     if floor(100*t)/100 >= t_record
%         % update marker positions in place (no delete/recreate)
%         set(m1, 'XData', x(1,:), 'YData', x(2,:), 'ZData', x(3,:)) ;
% 
%         % append this frame's position to each particle's trail
%         for k = 1:num
%             addpoints(trails(k), x(1,k), x(2,k), x(3,k)) ;
%         end
% 
%         set(titleHandle, 'String', ['t = ' num2str(round(t,1))]) ;
%         drawnow limitrate ;
% %         writeVideo(vidObj, getframe(figure(1)))
%         t_record = t_record + dt_record ;
%         count_record = count_record + 1 ;
%         % fileName = strcat(num2str(count_record), '.mat') ;
%         % save(fullfile(savdir, fileName), 'x', 'y', 't')
%     end   
%     for i = 1:num
%         for j = i+1:num
%             if norm(x(:,i) - x(:,j)) < 1e-4
%                 break ;
%             end
%         end
%     end
% end
% 
% close all
% % close(vidObj)
% 
% 
% %%
% 
% function v = velocity_field_nbody(x,i)
%     v = 0 ;
%     for j = 1:size(x,2)
%         if j ~= i
%             % Cr = 1.0;    lr = 0.5;   % repulsion strength & length scale (short range)
%             % Ca = 0.5;    la = 2.0;   % attraction strength & length scale (long range)
%             % r = norm(x(:,i) - x(:,j)) ;
%             % Wp = -Cr/lr*exp(-r/lr) + Ca/la*exp(-r/la) ;
%             % v = v - Wp/r * (x(:,i) - x(:,j)) ;
%             v = v + (x(:,i) - x(:,j))/ (norm(x(:,i) - x(:,j))^3) ;
%         end
%     end
% end