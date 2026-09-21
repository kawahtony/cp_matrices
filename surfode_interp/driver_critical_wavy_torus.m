clear ; clc ; close all ;

% Set random seed
rng(42)

%% Load surface files

load('data/wavy_torus_dx=0.025.mat') ;


%% Compute mean curvature (use unit outward normal)

% % Differentiation and interpolation matrices
[Dxc, Dyc, Dzc] = firstderiv_cen2_3d_matrices(x1d, y1d, z1d, band) ;
Emat = interp3_matrix(x1d, y1d, z1d, cpx_band, cpy_band, cpz_band, q, band);

H = wavyTorusMeanCurvature(uu_band, vv_band) ;

% H = Emat * H ; 


%% Define velocity field

Hx = Emat * (Dxc * H) ;
Hy = Emat * (Dyc * H) ; 
Hz = Emat * (Dzc * H) ;

kappa = 1 ; % ascent
% kappa = - 1 ;% descent

fx = kappa*Hx ; 
fy = kappa*Hy ;
fz = kappa*Hz ; 

% Project onto tangent space
for k = 1:length(band)
    ftemp = [fx(k); fy(k); fz(k)] ;
    ntemp = [nx(k); ny(k); nz(k)] ;
    ftemp = ftemp - dot(ftemp, ntemp)*ntemp ;
    fx(k) = ftemp(1) ; fy(k) = ftemp(2) ; fz(k) = ftemp(3) ;
end


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

% Number of points to be explored
numpt = 200 ;


% Get random indices of closest points
random_idx = randperm(length(band), numpt)' ;



%% Time-stepping

tol = 0.1 * dx^2 ;

critical_stored = zeros(3, numpt) ;
res_stored = zeros(1, numpt) ;
t_stored = zeros(1, numpt) ;


for k = 1:numpt

    idx = random_idx(k) ;
    x = [cpx_band(idx); cpy_band(idx); cpz_band(idx)] ;

    % Initialize CFL condition, time and residual
    cfl_constant = max(1, max(abs(fx)) + max(abs(fy)) + max(abs(fy))) ;
    t = 0 ; 
    count = 0 ;
    res = 1 ;
    Etemp = interp3_matrix(x1d, y1d, z1d, x(1), x(2), x(3), q, band) ;
    H_stored = Etemp*H ;
    avg_count = 10  ;

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

        % Update time and count 
        count = count + 1 ;
        t = t + dt ;

        % Update cfl constant
        cfl_constant = max(1, norm(Etemp*[fx, fy, fz], 1)) ; 


        % Check convergence  
        Etemp = interp3_matrix(x1d, y1d, z1d, x(1), x(2), x(3), q, band) ;
        Hval = Etemp*H ; 
        H_stored(count) = Hval ; 
        if count <= avg_count
            res = 1 ; 
        else    
            Havg0 = mean(H_stored(count - avg_count : count-1)) ;
            Havg = mean(H_stored(count - avg_count + 1 : count)) ;
            res = abs(Havg - Havg0)/max(abs(Havg), 1) ;
        end         


    end

    critical_stored(:,k) = x ;
    res_stored(k) = res ;
    t_stored(k) = t ; 
    disp([num2str(k), '/', num2str(numpt), ', t = ', num2str(t)])
    figure(1) ;
    plot3(x(1), x(2), x(3), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 7)

end



%% Save 

if kappa == 1
    save('critical_temp_ascent.mat', 'critical_stored', 'res_stored', 't_stored')
elseif kappa == -1
    save('critical_temp_descent.mat', 'critical_stored', 'res_stored', 't_stored')
end