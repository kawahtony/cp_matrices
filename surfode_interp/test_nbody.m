% %%  n-body problems on triangulated surfaces
% 
% clear ; close all ; clc ;
% 
% 
% % grid size & time-stepping method
% dx = 0.05 ;
% 
% % choice of triangulated surfaces
% PlyFile = 'bunny.ply' ;
% % shape = 'pig_loop2' ;
% % shape = 'cow' ;
% 
% switch PlyFile
%     case 'bunny.ply'
%         xmin = -1.5 ; xmax = 1.5 ;
%         ymin = -1.5 ; ymax = 1.5 ;
%         zmin = -1.5 ; zmax = 1.5 ;
%     case 'pig_loop2.ply'
%         xmin = -1.5 ; xmax = 1.5 ;
%         ymin = -1.5 ; ymax = 1.5 ;
%         zmin = -1.5 ; zmax = 1.5 ;
%     case 'cow.ply'
%         xmin = -1.5 ; xmax = 1.5 ;
%         ymin = -2.5 ; ymax = 2.5 ;
%         zmin = -0.5 ; zmax = 3.5 ;
% end
% 
% 
% ptCloud = pcread(PlyFile) ;
% vertices = double(ptCloud.Location) ;
% temp = vertices(:,2) ;
% vertices(:,2) = vertices(:,3) ; 
% vertices(:,3) = temp ; 
% clear temp ;
% 
% faces = MyRobustCrust(vertices) ;
% faces = double(faces) ;
% 
% 
% %%
% 
% 
% % make vectors of x, y, positions of the grid
% x1d = (xmin:dx:xmax)' ;
% y1d = (ymin:dx:ymax)' ;
% z1d = (zmin:dx:zmax)' ;
% 
% 
% % embedding space
% [xx, yy, zz] = meshgrid(x1d, y1d, z1d) ;
% 
% %%
% 
% dim = 3 ; % Dimension of embedding space
% q = 3 ; % Degree of interpolating polynomial
% 
% bw = 1 + sqrt(dim)*(q+1)/2 ;
% 
% [sdist, cpx, cpy, cpz, faceidx] = tri2sdf(xx, yy, zz, faces, vertices, bw*dx) ;
% 
% 
% %% Banding
% 
% cpx = cpx(:) ;
% cpy = cpy(:) ;
% cpz = cpz(:) ;
% 
% 
% band = find(abs(sdist) <= bw*dx) ;
% 
% % Store closest points in the band
% xx_band = xx(band) ;
% yy_band = yy(band) ;
% zz_band = zz(band) ;
% cpx_band = cpx(band) ;
% cpy_band = cpy(band) ;
% cpz_band = cpz(band) ;
% sdist_band = sdist(band) ;
% 
% 

% 
% save('data/bunny_temp.mat') ;


%%

clc ; clear ; close all ; 

jade = [0, 168, 107]/255 ; 

load('data/bunny_temp.mat') ;
% load('data/trefoil_knot_tube_dx=0.02.mat')

xx = xx(:) ; yy = yy(:) ; zz = zz(:) ;
xx_band = xx(band) ;
yy_band = yy(band) ;
zz_band = zz(band) ;


%% Compute unit outward normal

[Dx, Dy, Dz] = firstderiv_cen2_3d_matrices(x1d, y1d, z1d, band) ;

nx = (xx_band - cpx_band).*sign(sdist_band) ; 
ny = (yy_band - cpy_band).*sign(sdist_band)  ;
nz = (zz_band - cpz_band).*sign(sdist_band)  ; 
nnorm = sqrt(nx.^2 + ny.^2 + nz.^2) ;

phi_x = Dx*sdist_band ;
phi_y = Dy*sdist_band ;
phi_z = Dz*sdist_band ;
for k = 1:length(nnorm)
    if nnorm(k) < 1e-10
        nx(k) = phi_x(k) ; 
        ny(k) = phi_y(k) ;
        nz(k) = phi_z(k) ;
        nnorm(k) = sqrt(nx(k)^2 + ny(k)^2 + nz(k)^2) ;
    end
end

nx = nx./nnorm ; 
ny = ny./nnorm ; 
nz = nz./nnorm ;



%%

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

figure(1) ; clf ;
p = patch('Faces', faces, 'Vertices', vertices, 'FaceAlpha', 0.5, 'FaceColor', jade, ...
          'EdgeColor', 'none') ;

axis equal ; axis tight ; axis off ; hold on ;

% distinct color per particle, shared by its marker and its trail
colors = lines(num) ;

% one animatedline per particle for the trail history: addpoints()
% appends in place, so the trail never needs to be deleted/redrawn from
% scratch as it grows
trailMaxPts = 50 ;   % cap trail length so draw cost stays bounded over long runs
trails = gobjects(num,1) ;
for k = 1:num
    trails(k) = animatedline('Color', colors(k,:), 'LineWidth', 1, ...
        'MaximumNumPoints', trailMaxPts) ;
    addpoints(trails(k), x(1,k), x(2,k), x(3,k)) ;
end

% current particle positions; updated via set(...,'XData',...) each frame
% instead of delete()+replot
m1 = scatter3(x(1,:), x(2,:), x(3,:), 36, colors, 'filled', 'MarkerEdgeColor', 'k') ;

titleHandle = title('t = 0', 'FontSize', 16) ;

view(-181.5807, 17.7373)

k1_x = zeros(size(x)) ;
k1_v = zeros(size(x)) ;

alpha = 1 ; 
beta = 0.5 ;

while t < T

    fmax = 1 ;
    for k = 1:num
        fmax = max(fmax, norm(f(x,k))) ;
    end  

    dt = dx/fmax ; disp(dt)

    if T - t < dt
        dt = T - t ;
    end

    % stage 1

    for k = 1:num      
        Etemp = interp3_matrix(x1d, y1d, z1d, x(1,k), x(2,k), x(3,k), 3, band) ;
        ntemp = ( Etemp*[nx, ny, nz] )' ; 
        k1_x(:,k)= v(:,k) ;
        k1_v(:,k) = (alpha - beta*norm(v(:,k))^2)*v(:,k) + f(x,k) ;
        k1_v(:,k) = k1_v(:,k) - dot(k1_v(:,k), ntemp)*ntemp ;
        
    end
    % % stage 2
    % xtemp = x + dt*k1 ;
    % k2 = zeros(size(x)) ;
    % for k = 1:num              
    %     Etemp = interp3_matrix(x1d, y1d, z1d, xtemp(1,k), xtemp(2,k), xtemp(3,k), 3, band) ;
    %     ntemp = Etemp*[nx, ny, nz] ; ntemp = ntemp' ; 
    %     k2(:,k) = f(xtemp,k) ;
    %     k2(:,k) = k2(:,k) - dot(k2(:,k), ntemp)*ntemp ;             
    % end
    % 
    % % Update from RK2
    % x = x + dt*(k1/2 + k2/2) ;

    % Update from FE
    x = x + dt*k1_x ;
    v = v + dt*k1_v ;

    % Projection
    for k = 1:num
        Etemp =  interp3_matrix(x1d, y1d, z1d, x(1,k), x(2,k), x(3,k), 3, band) ;
        ntemp = ( Etemp*[nx, ny, nz] )' ;
        xproj = ( Etemp*[cpx_band, cpy_band, cpz_band] )' ;
        x(:,k) = xproj ;
        v(:,k) = v(:,k) - dot(v(:,k), ntemp)*ntemp ;
        % if norm(x(:,k) - xproj) > (0.1*bw)*dx
        %     x(:,k) = xproj ;
        % end
    end
    t = t + dt ;

    if floor(100*t)/100 >= t_record
        % update marker positions in place (no delete/recreate)
        set(m1, 'XData', x(1,:), 'YData', x(2,:), 'ZData', x(3,:)) ;

        % append this frame's position to each particle's trail
        for k = 1:num
            addpoints(trails(k), x(1,k), x(2,k), x(3,k)) ;
        end

        set(titleHandle, 'String', ['t = ' num2str(round(t,1))]) ;
        drawnow limitrate ;
%         writeVideo(vidObj, getframe(figure(1)))
        t_record = t_record + dt_record ;
        count_record = count_record + 1 ;
        % fileName = strcat(num2str(count_record), '.mat') ;
        % save(fullfile(savdir, fileName), 'x', 'y', 't')
    end   
    for i = 1:num
        for j = i+1:num
            if norm(x(:,i) - x(:,j)) < 1e-4
                break ;
            end
        end
    end
end

close all
% close(vidObj)


%%

function v = velocity_field_nbody(x,i)
    v = 0 ;
    for j = 1:size(x,2)
        if j ~= i
            % Cr = 1.0;    lr = 0.5;   % repulsion strength & length scale (short range)
            % Ca = 0.5;    la = 2.0;   % attraction strength & length scale (long range)
            % r = norm(x(:,i) - x(:,j)) ;
            % Wp = -Cr/lr*exp(-r/lr) + Ca/la*exp(-r/la) ;
            % v = v - Wp/r * (x(:,i) - x(:,j)) ;
            v = v + (x(:,i) - x(:,j))/ (norm(x(:,i) - x(:,j))^3) ;
        end
    end
end