%% "Gradient descent" to a desired curvature value

clear ; clc ; close all ;


%% Load surface files

load('data/wavy_torus_dx=0.025.mat') ;


%% Compute mean curvature (use unit outward normal)

% % Differentiation and interpolation matrices
% Lmat = laplacian_3d_matrix(x1d, y1d, z1d, 2, band);
[Dxc, Dyc, Dzc] = firstderiv_cen2_3d_matrices(x1d, y1d, z1d, band) ;
Emat = interp3_matrix(x1d, y1d, z1d, cpx_band, cpy_band, cpz_band, q, band);
% 
% % Discretized laplacian of closest point functions
% Lcpx = Emat * (Lmat * cpx_band) ; 
% Lcpy = Emat * (Lmat * cpy_band) ; 
% Lcpz = Emat * (Lmat * cpz_band) ; 
% 
% H = zeros(size(band)) ;
% for k = 1:length(band)
%     H(k) = -dot([Lcpx(k); Lcpy(k); Lcpz(k)], [nx(k); ny(k); nz(k)]) ;
% end
% 
% H = Emat * H ; 

H = wavyTorusMeanCurvature(uu_band, vv_band) ;
H = Emat * H ; 

%% Define velocity field

Hx = Emat * (Dxc * H) ;
Hy = Emat * (Dyc * H) ; 
Hz = Emat * (Dzc * H) ;

kappa = 1 ; 

fx = kappa*Hx ; 
fy = kappa*Hy ;
fz = kappa*Hz ; 

% % Project onto tangent space
% for k = 1:length(band)
%     ftemp = [fx(k); fy(k); fz(k)] ;
%     ntemp = [nx(k); ny(k); nz(k)] ;
%     ftemp = ftemp - dot(ftemp, ntemp)*ntemp ;
%     fx(k) = ftemp(1) ; fy(k) = ftemp(2) ; fz(k) = ftemp(3) ;
% end


%% Initialize plot 

[up, vp] = meshgrid(linspace(0,2*pi, 300)) ;

xp  = ((169*cos(2*up + 3*vp)/2500 + 13/25).*cos(vp) + 6*cos(5*up)/25 + 41/20).*cos(up) ;
          
yp = ((169*cos(2*up + 3*vp)/2500 + 13/25).*cos(vp) + 6*cos(5*up)/25 + 41/20).*sin(up) ;

zp = (169*cos(2*up + 3*vp)/2500 + 13/25).*sin(vp) + 9*sin(5*up)/50 ;


Eplot = interp3_matrix(x1d, y1d, z1d, xp(:), yp(:), zp(:), q, band) ;

figure(1) ; clf ; hold on ;

Hp = Eplot*H ;
Hp = reshape(Hp, size(xp)) ;
surf(xp, yp, zp, Hp , 'FaceAlpha', 0.9) ;
shading interp
axis equal 
view(3)



%% Time discretization

% Simulation time
T = 20 ;

% Initial location
% x = [1.2721; -1.6560; -0.4047] ; % descent
% x = [1.83752; 0; -0.0142673] ; % ascent
% x = [-2.10874; -0.321765; 0.232841] ; % ascent
x = [0.7123; 2.1638; 0.5530] ;

figure(1) ; 
plot3(x(1), x(2), x(3), 'ro', 'MarkerFaceColor', 'r')

x_stored = [] ; x_stored(:,1) = x ;
t_stored = 0 ;

Etemp = interp3_matrix(x1d, y1d, z1d, x(1), x(2), x(3), q, band) ;

H_stored = Etemp*H ;
Havg_stored = Etemp*H ;
avg_count = 20  ;

count = 0 ; 

% Initialize CFL condition, time and residual
cfl_constant = max(1, max(abs(fx)) + max(abs(fy)) + max(abs(fy))) ;
t = 0 ;
res = 1 ; 
tol = 0.01*dx^2 ;


%% Time-stepping


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

    %% Update time and save
    t = t + dt ;
    count = count + 1 ; 
    x_stored(:,count) = x ; 
    t_stored(count) = t ; 

    % Check convergence
    if count <= avg_count
        res = 1 ; 
    else
        Etemp = interp3_matrix(x1d, y1d, z1d, x(1), x(2), x(3), q, band) ;
        Hval = Etemp*H ; 
        H_stored(count) = Hval ; 
        Havg0 = mean(H_stored(count - avg_count : count-1)) ;
        Havg = mean(H_stored(count - avg_count + 1 : count)) ;
        Havg_stored(count) = Havg ; 
        res = abs(Havg - Havg0)/max(abs(Havg), 1) ; disp(res)
    end    

    disp(['t = ', num2str(t)])

end


figure(1) ;
plot3(x_stored(1,:), x_stored(2,:), x_stored(3,:), 'r-', 'LineWidth', 1.5)


figure(2) ; 
plot3(x_stored(1,1), x_stored(2,1), x_stored(3,1), 'ro', 'MarkerFaceColor', 'r') ; hold on ;
plot3(x_stored(1,:), x_stored(2,:), x_stored(3,:), 'r-', 'LineWidth', 1.5)


figure(3) ; 
plot(t_stored, H_stored, 'r-', 'LineWidth', 1.5)


figure(4) ; 
plot(t_stored, Havg_stored, 'r-', 'LineWidth', 1.5)