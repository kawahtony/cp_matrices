%% Test code to find the "period" of contour

clear ; clc ; close all ;


%% Load surface data

surface_data = 'data/twisted_elliptic_torus_dx=0.025.mat' ;
load(surface_data) ;


%% Compute mean curvature (use unit outward normal)

% Differentiation and interpolation matrices
Lmat = laplacian_3d_matrix(x1d, y1d, z1d, 2, band);
[Dxc, Dyc, Dzc] = firstderiv_cen2_3d_matrices(x1d, y1d, z1d, band) ;
Emat = interp3_matrix(x1d, y1d, z1d, cpx_band, cpy_band, cpz_band, q, band);

% Discretized laplacian on closest point functions
Lcpx = Emat * (Lmat * cpx_band) ; 
Lcpy = Emat * (Lmat * cpy_band) ; 
Lcpz = Emat * (Lmat * cpz_band) ; 


H = zeros(size(band)) ;
for k = 1:length(band)
    H(k) = -dot([Lcpx(k); Lcpy(k); Lcpz(k)], [nx(k); ny(k); nz(k)]) ;
end

H = Emat * H ; 


%% Define velocity field

Hx = Emat * (Dxc * H) ;
Hy = Emat * (Dyc * H) ; 
Hz = Emat * (Dzc * H) ;

fx = zeros(size(Hx)) ; fy = fx ; fz = fx ; 

% Velocity along contour
for k = 1:length(fx)
    vtemp = [Hx(k); Hy(k); Hz(k)] ;
    ntemp = [nx(k); ny(k); nz(k)] ;
    vtemp = cross(vtemp, ntemp) ; vtemp = vtemp / norm(vtemp) ;
    fx(k) = vtemp(1) ;
    fy(k) = vtemp(2) ;
    fz(k) = vtemp(3) ;
end


%% Initialize plot 

xp = vertices(:,1) ;
yp = vertices(:,2) ;
zp = vertices(:,3) ;

Eplot = interp3_matrix(x1d, y1d, z1d, xp, yp, zp, q, band) ;

figure(1) ; clf ; hold on ;
trisurf(faces, xp, yp, zp, Eplot*H, 'FaceAlpha', 0.5) ;
shading interp
axis equal 
view(3)


%% Time discretization

% Initialize time 
t = 0 ;

% Initialize CFL condition
cfl_constant = max(1, max(abs(fx)) + max(abs(fy)) + max(abs(fy))) ;


contour_id = 3 ; 

switch contour_id
    case 1
        x0 = [1.52133; 0.96879; 0.121337] ;
    case 2
        x0 = [1.43633; 1.05154; 0.326361] ;
    case 3
        x0 = [1.45038; 1.06468; 0.411222] ;
    case 4
        x0 = [1.15183; 1.17467; 0.307772] ;
end


figure(1) ; plot3(x0(1), x0(2), x0(3), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 7)


%% Time-stepping

x = x0  ;
x_stored = [] ; t_stored = [] ; H_stored = [] ;
x_stored(:,1) = x ; t_stored(1) = 0 ; 
Etemp = interp3_matrix(x1d, y1d, z1d, x(1), x(2), x(3), q, band) ; H_stored(1) = Etemp*H ;
count = 1 ;
res = 1 ;

while count <= 20 || res > 0.5*dx 

    dt = dx / cfl_constant ;

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
    x = xproj ;

    % Update cfl constant
    cfl_constant = max(1, norm(Etemp*[fx, fy, fz], 1)) ; 

    % Save the second position
    if count == 1
        x1 = x ;
    end

    % Update residual to the initial position
    res = min(norm(x - x0), norm(x - x1)) ;
    Etemp = interp3_matrix(x1d, y1d, z1d, x(1), x(2), x(3), q, band) ;
    Htemp = Etemp*H ;
    disp([t, res, Htemp])

    %% Save solution
    count = count + 1 ;
    t = t + dt ;
    x_stored(:,count) = x ;  
    t_stored(count) = t ;
    H_stored(count) = Htemp ;

end


%% Plot  

figure(1) ; 
plot3(x_stored(1,:), x_stored(2,:), x_stored(3,:), 'r-', 'LineWidth', 1.5)
axis equal 

figure(2) ; clf ; 
plot3(x_stored(1,:), x_stored(2,:), x_stored(3,:), 'r-', 'LineWidth', 1.5)


figure(3) ; 
plot(t_stored, H_stored, 'k-', 'LineWidth', 1.5)


%% Save

save(['twisted_torus_contour_', num2str(contour_id), '.mat'], 'x_stored', 't_stored')