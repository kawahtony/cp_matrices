%% n-body on ellipsoid

clc ; clear ; close all ;

azure = [0, 128, 255]/255 ;
coral = [255, 127, 80]/255 ;
jade = [0, 168, 107]/255 ;

format long


%% Load closest point of trefoil knot tubulbar surface

a = 1.5; % length along z
b = 1; % length along xy plane

dx = 0.08 ;

x1d = (-2.01:dx:2.01)';
y1d = x1d;
z1d = x1d;

[xx,yy,zz] = meshgrid(x1d, y1d, z1d);

[cpx,cpy,cpz,dist] = cpEllipsoid(xx, yy, zz, [a b], [], 'z');


%% Banding

dim = 3 ;
q = 3 ;
order = 2 ; 
bw = 1.0001*sqrt((dim-1)*((q+1)/2)^2 + ((order/2+(q+1)/2)^2));

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

%% Differentiation matrices

[Dxc, Dyc, Dzc] = firstderiv_cen2_3d_matrices(x1d, y1d, z1d, band) ;
Emat = interp3_matrix(x1d, y1d, z1d, cpx_band, cpy_band, cpz_band, q, band);


%% Define velocity field

xs = -0.38897 ; 
ys = -0.767841 ;
zs = 0.763562 ;

temp = -((xx_band - xs).^2 + (yy_band - ys).^2 + (zz_band - zs).^2) ;
g = exp(temp/0.5^2) ;

g = Emat * g ;

gx = Emat * (Dxc * g) ;
gy = Emat * (Dyc * g) ;
gz = Emat * (Dzc * g) ;

kappa = 1 ; 

fx = kappa * gx ; 
fy = kappa * gy ; 
fz = kappa * gz ; 

% % Project onto tangent space
% for k = 1:length(band)
%     ftemp = [fx(k); fy(k); fz(k)] ;
%     ntemp = [nx(k); ny(k); nz(k)] ;
%     ftemp = ftemp - dot(ftemp, ntemp)*ntemp ;
%     fx(k) = ftemp(1) ; fy(k) = ftemp(2) ; fz(k) = ftemp(3) ;
% end


%% Initialize plot

[xp, yp, zp] = paramEllipsoid(100, [a b], [], 'z');

Eplot = interp3_matrix(x1d, y1d, z1d, xp(:), yp(:), zp(:), q, band) ;


figure(1) ; clf ;
surf(xp, yp, zp, reshape(Eplot*g, size(xp)), 'EdgeColor', 'None', 'FaceAlpha', 0.6); hold on ;
colormap('parula')
shading interp
axis equal;


%% Time discretization

% Simulation time
T = 10 ;

x = [-0.753739; -0.237113; 0.919361] ;

plot3(x(1), x(2), x(3), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r')


% Initialize CFL condition, time and residual
cfl_constant = max(1, max(abs(fx)) + max(abs(fy)) + max(abs(fy))) ;
t = 0 ;
res = 1 ; 
tol = 0.01*dx^2 ;

count = 0 ;

x_stored = x ; 
t_stored = 0 ;

Etemp = interp3_matrix(x1d, y1d, z1d, x(1), x(2), x(3), q, band) ;
g_stored = Etemp * g ;
gavg_stored = Etemp*g ; 
avg_count = 10 ; 



%% Time stepping

while t < T && res > tol 

    x0 = x ;
    dt = dx / cfl_constant ;

    if T - t < dt
            dt = T - t ;
    end
    
    %% Explicit Trapezoidal (RK2) time integration

    % Stage 1
    Etemp = interp3_matrix(x1d, y1d, z1d, x(1), x(2), x(3), q, band) ;
    k1 = ( Etemp*[fx, fy, fz] )' ; 

    % Stage 2
    xtemp = x + dt*k1 ;
    Etemp = interp3_matrix(x1d, y1d, z1d, xtemp(1), xtemp(2), xtemp(3), q, band) ; 
    k2 = ( Etemp*[fx, fy, fz] )' ; 

    % Update from TVD-RK2
    x = x + dt*(k1/2 + k2/2) ;        

    %% Final Projection
    Etemp = interp3_matrix(x1d, y1d, z1d, x(1), x(2), x(3), q, band) ;
    xproj = ( Etemp*[cpx_band, cpy_band, cpz_band] )' ; 
    if norm(x - xproj) > dx/sqrt(2)
        x = xproj ; 
    end

    % Update cfl constant
    cfl_constant = max(1, norm(Etemp*[fx, fy, fz], 1)) ; 

    %% Update time and count
    t = t + dt ;
    count = count + 1 ;
    x_stored(:,count) = x ; 
    t_stored(count) = t ; 


    % Check convergence
    if count <= avg_count
        res = 1 ; 
    else
        Etemp = interp3_matrix(x1d, y1d, z1d, x(1), x(2), x(3), q, band) ;
        gval = Etemp*g ; 
        g_stored(count) = gval ; 
        gavg0 = mean(g_stored(count - avg_count : count-1)) ;
        gavg = mean(g_stored(count - avg_count + 1 : count)) ;
        gavg_stored(count) = gavg ; 
        res = abs(gavg - gavg0)/max(abs(gavg), 1) ; disp(res)
    end    


end

plot3(x_stored(1,:), x_stored(2,:), x_stored(3,:), 'r-', 'LineWidth', 1.5)

diff = norm([x_stored(1,end) - xs, x_stored(2,end) - ys, x_stored(3,end) - zs], 2) ;
disp(diff)
