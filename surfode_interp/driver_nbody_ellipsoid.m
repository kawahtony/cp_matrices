%% n-body on ellipsoid

clc ; clear ; close all ;

rng(42)

azure = [0, 128, 255]/255 ;
coral = [255, 127, 80]/255 ;
jade = [0, 168, 107]/255 ;


%% Load closest point of trefoil knot tubulbar surface

a = 1.5; % length along z
b = 1; % length along xy plane

dx = 0.05 ;

x1d = (-2.01:dx:2.01)';
y1d = x1d;
z1d = x1d;

[xx,yy,zz] = meshgrid(x1d, y1d, z1d);

[cpx,cpy,cpz,dist] = cpEllipsoid(xx, yy, zz, [a b], [], 'z');


%% Banding

dim = 3 ;
q = 3 ;
bw = 1 + sqrt(dim)*(q + 1)/2 ;

cpx = cpx(:) ; cpy = cpy(:) ; cpz = cpz(:) ;
xx = xx(:) ; yy = yy(:) ; zz = zz(:) ;
dist = dist(:) ;

band = find(abs(dist) <= bw*dx);

xx_band = xx(band) ;
yy_band = yy(band) ;
zz_band = zz(band);
cpx_band = cpx(band) ;
cpy_band = cpy(band) ;
cpz_band = cpz(band) ;


%% Compute unit outward normal

u2 = acos(cpz_band/a) ;
th2 = atan2(cpx_band, cpy_band) ;

nx = zeros(size(band)) ; ny = nx ; nz = nx ;
for k = 1:length(nx)
    r_u2 = [b*cos(u2(k))*sin(th2(k)); b*cos(u2(k))*cos(th2(k)); -a*sin(u2(k))] ;
    r_th = [b*sin(u2(k))*cos(th2(k)); -b*sin(u2(k))*sin(th2(k)); 0] ;
    ntemp = cross(r_u2, r_th) ; ntemp = ntemp / norm(ntemp) ;
    nx(k) = ntemp(1) ; ny(k) = ntemp(2) ; nz(k) = ntemp(3) ;
end


%%

numpt = 10 ;
th0 = linspace(0, 2*pi*(1 - 1/numpt), numpt)' ;
phi0 = pi/5 ;

x = [b*sin(phi0)*sin(th0), b*sin(phi0)*cos(th0), a*cos(phi0)*ones(numpt,1)]' ;

%% plot

[xp, yp, zp] = paramEllipsoid(100, [a b], [], 'z');

figure(1) ; clf ;
surf(xp,yp,zp,zp, 'EdgeColor', 'None', 'FaceAlpha', 0.5); hold on ;
colormap('parula')
shading interp
axis equal;


% distinct color per particle, shared by its marker and its trail
colors = lines(numpt) ;

m1 = scatter3(x(1,:), x(2,:), x(3,:), 36, colors, 'filled', 'MarkerEdgeColor', 'k') ;
titleHandle = title('t = 0', 'FontSize', 16) ;



%% Set-up for n-body simulation

% define velocity field
f = @velocity_field_nbody ;


% time-stepping initialization
t = 0 ;
T = 100 ;
v = randn(3,numpt) ;


dt_record = 0.01 ;
t_record = dt_record ;


%% Initialize plot


k1_x = zeros(size(x)) ;
k1_v = zeros(size(x)) ;

% sliding window of the most recent trail points (fixed length, oldest
% points drop off as new ones are recorded)
trailLen = 50 ;
trail_x = nan(numpt, trailLen) ;
trail_y = nan(numpt, trailLen) ;
trail_z = nan(numpt, trailLen) ;

trail_x(:,end) = x(1,:)' ;
trail_y(:,end) = x(2,:)' ;
trail_z(:,end) = x(3,:)' ;

% one trail line per particle, colored to match its marker
trailHandles = gobjects(numpt,1) ;
for k = 1:numpt
    trailHandles(k) = plot3(trail_x(k,:), trail_y(k,:), trail_z(k,:), ...
        '-', 'Color', colors(k,:), 'LineWidth', 1) ;
end

while t < T

    fmax = 1 ;
    for k = 1:numpt
        fmax = max(fmax, norm(f(x,k))) ;
    end

    dt = dx/fmax ;

    if T - t < dt
        dt = T - t ;
    end

    for k = 1:numpt
        Etemp = interp3_matrix(x1d, y1d, z1d, x(1,k), x(2,k), x(3,k), q, band) ;
        ntemp = ( Etemp*[nx, ny, nz] )' ;
        k1_x(:,k)= v(:,k) ;
        k1_v(:,k) = f(x,k) ;
        k1_v(:,k) = k1_v(:,k) - dot(k1_v(:,k), ntemp)*ntemp ;
    end

    % Update from FE
    x = x + dt*k1_x ;
    v = v + dt*k1_v ;

    % Projection
    for k = 1:numpt
        Etemp =  interp3_matrix(x1d, y1d, z1d, x(1,k), x(2,k), x(3,k), q, band) ;
        ntemp = ( Etemp*[nx, ny, nz] )' ;
        xproj = ( Etemp*[cpx_band, cpy_band, cpz_band] )' ;
        x(:,k) = xproj ;
        v(:,k) = v(:,k) - dot(v(:,k), ntemp)*ntemp ;
    end
    t = t + dt ;

    if t >= t_record
        % update marker positions in place (no delete/recreate)
        set(m1, 'XData', x(1,:), 'YData', x(2,:), 'ZData', x(3,:)) ;

        % slide the trail window forward and redraw for each particle
        trail_x = [trail_x(:,2:end), x(1,:)'] ;
        trail_y = [trail_y(:,2:end), x(2,:)'] ;
        trail_z = [trail_z(:,2:end), x(3,:)'] ;
        for k = 1:numpt
            set(trailHandles(k), ...
                'XData', trail_x(k,:), ...
                'YData', trail_y(k,:), ...
                'ZData', trail_z(k,:)) ;
        end

        drawnow limitrate ;
        t_record = t_record + dt_record ;
        set(titleHandle, 'String', ['t = ' num2str(round(t,1))]) ;
    end

    for i = 1:numpt
        for j = i+1:numpt
            if norm(x(:,i) - x(:,j)) < 1e-4
                break ;
            end
        end
    end
end


% x0 = [-0.1709    0.7277   -0.4183   -0.0991;
%     0.6306   -0.3386   -0.4685   -0.0892;
%     1.1355    0.8948    1.1672    1.4866] ;
%
% figure(2) ; clf ;
% surf(xp,yp,zp,zp, 'EdgeColor', 'None', 'FaceAlpha', 0.5); hold on ;
% plot3(x0(1,:), x0(2,:), x0(3,:), 'ro', 'MarkerFaceColor', 'r') ;
% shading interp
% axis equal;


%%

function v = velocity_field_nbody(x,i)
    v = 0 ;
    for j = 1:size(x,2)
        if j ~= i
            v = v + (x(:,i) - x(:,j))/ (norm(x(:,i) - x(:,j))^3) ;
        end
    end
end