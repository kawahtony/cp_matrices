clear ; clc ; close all ;


%% Discretize embedding space and compute closest point of square

dx = 0.05 ; 

x1d = (-1.5:dx:1.5)' ; 
y1d = x1d ;

[xx, yy] = meshgrid(x1d, y1d) ;

[cpx, cpy, sdist] = cpSquare(xx, yy) ;


%% Banding 

dim = 2 ;
q = 3 ;
bw = 1 + sqrt(dim)*(q + 1)/2 ;


cpx = cpx(:) ; cpy = cpy(:) ;
xx = xx(:) ; yy = yy(:) ;
sdist = sdist(:) ;

band = find(abs(sdist) <= bw*dx) ;

xx_band = xx(band) ; 
yy_band = yy(band) ;
cpx_band = cpx(band) ; 
cpy_band = cpy(band) ;
sdist_band = sdist(band) ;


%% Compute unit outward normal

[Dxc, Dyc] = firstderiv_cen2_2d_matrices(x1d, y1d, band) ;


nx = (xx_band - cpx_band) .* sign(sdist_band) ;
ny = (yy_band - cpy_band) .* sign(sdist_band) ;

nnorm = sqrt(nx.^2 + ny.^2) ;

phi_x = Dxc * sdist_band ;
phi_y = Dyc * sdist_band ;
for k = 1:length(nnorm)
    if nnorm(k) < 1e-10
        nx(k) = phi_x(k) ; 
        ny(k) = phi_y(k) ;
        nnorm(k) = sqrt(nx(k)^2 + ny(k)^2) ;
    end
end

nx = nx./nnorm ; 
ny = ny./nnorm ; 


%% Velocity field

f1 = -cpy_band ; 
f2 = cpx_band ;


% Projection on tangent space
for k = 1:length(band)
    vtemp = [f1(k); f2(k)] ;
    ntemp = [nx(k); ny(k)] ;
    vtemp = vtemp - dot(vtemp, ntemp)*ntemp ;
    f1(k) = vtemp(1) ; f2(k) = vtemp(2) ;
end


%% Time discretization

% Initialize time 
t = 0 ;

% Initialize CFL condition
cfl_constant = max(abs(f1)) + max(abs(f2)) ;

% Simulation time
T = 2 ;

% Initial location
x = [0; 1] ;


%% Initialize plot


figure(1) ; clf ; hold on ;
plot([-1 1], [1 1], 'k:', 'LineWidth', 1.5)
plot([1 1], [-1 1], 'k:', 'LineWidth', 1.5)
plot([-1 1], [-1 -1], 'k:', 'LineWidth', 1.5)
plot([-1 -1], [-1 1], 'k:', 'LineWidth', 1.5)
plot(x(1), x(2), 'ro', 'MarkerFaceColor', 'r')
axis equal 
axis([-1.2 1.2 -1.2 1.2])



%% Time stepping

k1 = zeros(2,1) ; 
k2 = zeros(2,1) ;
count_proj = 0 ; 
count = 1 ;
x_stored = x' ;


while t < T

    dt = dx / cfl_constant ;

    if T - t < dt 
        dt = T - t ;  
    end

    %% Explicit Trapezoidal (RK2) time integrator

    % Stage 1 
    Etemp = interp2_matrix(x1d, y1d, x(1), x(2), q, band) ;
    k1 = [Etemp*f1 ; Etemp*f2] ;

    % Stage 2
    xtemp = x + dt*k1 ;  
    Etemp = interp2_matrix(x1d, y1d, xtemp(1), xtemp(2), q, band) ;
    k2 = [Etemp*f1 ; Etemp*f2] ;

    % Update solution
    x = x + dt*(k1/2 + k2/2) ;


    % Projection
    Etemp = interp2_matrix(x1d, y1d, x(1), x(2), q, band) ;
    xproj = [Etemp*cpx_band ; Etemp*cpy_band] ;
    if norm(x - xproj) > dx/sqrt(2)
        x = xproj ; 
        count_proj = count_proj + 1 ;
    end

    % Update cfl constant 
    cfl_constant = abs(Etemp*f1) + abs(Etemp*f2) ;

    % Update number of time steps
    count = count + 1 ;

    % Update time
    t = t + dt ;

    % Save numerical solution
    x_stored(count, :) = x' ;


end

plot(x_stored(:,1), x_stored(:,2), 'r-', 'LineWidth', 2)