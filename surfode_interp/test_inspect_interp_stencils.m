%% Heat equation on a circle

% This example solves the heat equation on the surface of a sphere,
% with initial conditions u = cos(4*theta), using *penalty approach*.

jade = [0, 168, 107]/255 ;
coral = [255, 127, 80]/255 ;

%% Construct a grid in the embedding space

dx = 0.02 ;

% make vectors of x,y, positions of the grid
x1d = (-2:dx:2)' ;
y1d = x1d ;

nx = length(x1d) ;
ny = length(y1d) ;


%% Find closest points on the surface
% For each point (x,y), we store the closest point on the circle
% (cpx,cpy)

% meshgrid is only needed for finding the closest points, not afterwards
[xx, yy] = meshgrid(x1d, y1d) ;
% function cpCircle for finding the closest points on a circle
[cpx, cpy, dist] = cpCircle(xx, yy) ;
% make into vectors
cpxg = cpx(:) ;
cpyg = cpy(:) ;


%% Banding : do calculations in a narrow band around the circle
dim = 2 ;    % dimension
p = 1 ;      % interpolation degree of differential operator
q = 3 ;      % interpolation degree of penalty term
order = 2 ;  % Laplacian order
% "band" is a vector of the indices of the points in the computation band.
% The formula for bw is found in [Ruuth & Merriman 2008] and the 1.0001 is
% a safety factor.
% bw = 1.0001*sqrt((dim-1)*((max(p,q)+1)/2)^2 + ((order/2+(max(p,q)+1)/2)^2)) ;
bw = (1 + sqrt(2*((q + 1)/2)^2)) ;
band = find(abs(dist) <= bw*dx) ;

% store closest points in the band
cpxg = cpxg(band) ;
cpyg = cpyg(band) ;
dist = dist(band) ;
xg = xx(band) ;
yg = yy(band) ;


% xi = 0.80778 ;
% yi = 0.589484 ; 
xi = 0.92 ;
yi = 0.44 ;

[Ibpt, Xgrid] = findGridInterpBasePt_vec({xi yi}, q, [x1d(1) y1d(1)], [dx dx]);


figure(1) ; clf ; 
plot(exp(1i*linspace(0,2*pi,300)), 'r-') ; hold on 
plot((1 + dx)*exp(1i*linspace(0,2*pi,300)), 'r-') ;
plot((1 - dx)*exp(1i*linspace(0,2*pi,300)), 'r-') ;
plot((1 + bw*dx)*exp(1i*linspace(0,2*pi,300)), '-', 'Color', jade) ;
plot((1 - bw*dx)*exp(1i*linspace(0,2*pi,300)), '-', 'Color', jade) ;
plot(xg, yg, '.', 'Color', jade)
plot(xi, yi, 'r*')
for i = 0:3
    for j = 0:3
        plot(Xgrid{1} + i*dx, Xgrid{2} + j*dx, 'ro', 'MarkerSize', 10)
    end
end
axis equal ;


% %% Function u in the embedding space
% % u is a function defined on the grid (e.g. temperature if solving heat
% % equation)
% 
% % assign some initial value (using initial value of cos (8*theta))
% [th,~] = cart2pol(xx,yy) ;
% u = sin(2*th) ;
% 
% % this makes u into a vector, containing only points in the band
% u = u(band) ;
% initialu = u ;       % store initial value
% 
% 
% %% Construct an interpolation matrix for closest point
% % This creates a matrix which interpolates data from the grid x1d, y1d,
% % onto the points [cpx, cpy]
% 
% disp('Constructing interpolation and laplacian matrices') ;
% 
% Ep = interp2_matrix(x1d, y1d, cpxg, cpyg, p, band) ;
% Eq = interp2_matrix(x1d, y1d, cpxg, cpyg, q, band) ;
% 
% % e.g. closest point extension : 
% % u = Ep*u (order p)
% % u = Eq*u (order q)
% 
% 
% %% Create Laplacian matrix
% 
% L = laplacian_2d_matrix(x1d, y1d, order, band) ;
% 
% 
% %% Construct an interpolation matrix for plotting on unit circle
% 
% % plotting grid on circle, using theta as a parametrization
% thetas = linspace(0,2*pi,100)' ;
% r = ones(size(thetas)) ;
% % plotting grid in Cartesian coords
% [xp, yp] = pol2cart(thetas,r) ;
% xp = xp(:) ;
% yp = yp(:) ;
% Eplot = interp2_matrix(x1d, y1d, xp, yp, p, band) ;
% 
% 
% %% Time-stepping for the heat equation
% 
% 
% Tf = 0.5 ;
% gamma = 2 * dim / dx^2 ;
% I = eye(size(L)) ;
% M = Ep*L - gamma*(I - Eq) ;
% 
% dt = 1/4*dx^2 ;
% 
% numtimesteps = ceil(Tf/dt) ;
% % adjust for integer number of steps
% dt = Tf/numtimesteps ;
% 
% 
% for kt = 1 : numtimesteps
% 
%     u = u + dt*M*u ;       
%     t = kt*dt ;
% 
%     if ( (mod(kt,100) == 0) || (kt < 10) || (kt == numtimesteps) )
%         figure(1) ; clf ;
%         plot2d_compdomain(u, xg, yg, dx, dx, 1) ;
%         title('Numerical soln in embedded domain') ;
%         xlabel('x') ;
%         ylabel('y') ;
%         plot(xp, yp, 'k-', 'linewidth', 2) ;
%     end
% end