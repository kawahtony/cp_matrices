% Euler's equations on sphere

clear
close all

% grid size
dx = 0.05 ;

% interpolation order
q = 3 ;

% make vectors of x, y, positions of the grid
x1d = (-2:dx:2)' ;
y1d = x1d ;
z1d = x1d ;


% embedding space
[xx, yy, zz] = meshgrid(x1d, y1d, z1d) ;
[cpx, cpy, cpz, sdist] = cpSphere(xx, yy, zz) ;
[face,vertex] = isosurface(xx, yy, zz, sdist, 0) ;
cpxg = cpx(:) ;
cpyg = cpy(:) ;
cpzg = cpz(:) ;
% sdist = sdist(:) ;


% band
dim = 3 ;
bw = 2*sqrt(dim) ;
band = find(abs(sdist) <= bw*dx) ;
band_ext = find(abs(sdist) <= (bw+2*sqrt(dim))*dx) ; % need an extended band to compute normals


% store closest points in the band
sdist_band = sdist(band) ;
sdist_band_ext = sdist(band_ext) ;
xx_band = xx(band) ;
yy_band = yy(band) ;
zz_band = zz(band) ;
cpx_band = cpxg(band) ;
cpy_band = cpyg(band) ;
cpz_band = cpzg(band) ;


% discrete cp-extension
E = interp3_matrix(x1d, y1d, z1d, cpx_band, cpy_band, cpz_band, q, band) ;


% compute outward unit normals
% [nx,ny,nz] = normals_from_cp(xg_band, yg_band, zg_band, ...
%     cpx_band, cpy_band, cpz_band, sdist_band, dx) ;
[Dxc, Dyc, Dzc] = firstderiv_cen2_3d_matrices(x1d, y1d, z1d, band, band_ext) ;
nx = E*( Dxc*sdist_band_ext ) ;
ny = E*( Dyc*sdist_band_ext ) ;
nz = E*( Dzc*sdist_band_ext ) ;

% npts = 40 ; normal = [1;2;3] ;
% PT = great_circle_points(npts, normal);
% idx =(floor(npts/3):npts) ;
% PT = PT(idx,:) ;
t = linspace(pi/3, 2*pi/3, 40)' ;
PT = [cos(t).*cos(t), cos(t).*sin(t), sin(t)] ;


%% Check residual of geodesic BVP at sampling points on great circle arc

% figure(1) ; clf ;  hold on
% patch('Faces', face, 'Vertices', vertex, 'FaceVertexCData', ...
%     0.8*ones(size(face)), 'EdgeColor', 'none', 'FaceColor','flat') ;
% plot3(PT(:,1), PT(:,2), PT(:,3), 'r-o', 'MarkerSize', 3, 'MarkerFaceColor', 'r')
% axis equal ;


% % Get projection matrices at sampling points
% 
% P11 = zeros(size(xx_band)) ;
% P12 = zeros(size(xx_band)) ;
% P13 = zeros(size(xx_band)) ;
% P21 = zeros(size(xx_band)) ;
% P22 = zeros(size(xx_band)) ;
% P23 = zeros(size(xx_band)) ;
% P31 = zeros(size(xx_band)) ;
% P32 = zeros(size(xx_band)) ;
% P33 = zeros(size(xx_band)) ;
% 
% for k = 1:length(xx_band)
%     ntemp = [nx(k); ny(k); nz(k)] ;
%     temp = eye(3) - ntemp*ntemp' ;
%     P11(k) = temp(1,1) ;
%     P12(k) = temp(1,2) ;
%     P13(k) = temp(1,3) ;
%     P21(k) = temp(2,1) ;
%     P22(k) = temp(2,2) ;
%     P23(k) = temp(2,3) ;
%     P31(k) = temp(3,1) ;
%     P32(k) = temp(3,2) ; 
%     P33(k) = temp(3,3) ;
% end


% res = zeros(size(PT,1)-2, 3) ;
% resnorm = zeros(size(PT,1)-2, 1) ;
% for k = 2:size(PT,1)-1
%     Etemp = interp3_matrix(x1d, y1d, z1d, PT(k,1), PT(k,2), PT(k,3), 3, band) ;
%     Ptemp = [Etemp*P11 Etemp*P12 Etemp*P13 ;
%              Etemp*P21 Etemp*P22 Etemp*P23 ;
%              Etemp*P31 Etemp*P32 Etemp*P33] ;    
%     res(k,:) = Ptemp*(PT(k+1,:) - 2*PT(k,:) + PT(k-1,:))' ;
%     resnorm(k) = norm(res(k)) ;
% end


%% Test direct minimization of arclength energy 


% PT(2:end-1,:) = PT(2:end-1,:) + 0.01*randn(size(PT(2:end-1,:))) ;
% 
% for k = 1:size(PT,1)
%     PT(k,:) = PT(k,:)/norm(PT(k,:)) ;
% end
    
% 
% figure(1) ; clf ;  hold on
% patch('Faces', face, 'Vertices', vertex, 'FaceVertexCData', ...
%     0.8*ones(size(face)), 'EdgeColor', 'none', 'FaceColor','flat') ;
% p = plot3(PT(:,1), PT(:,2), PT(:,3), 'r-o', 'MarkerSize', 3, 'MarkerFaceColor', 'r') ;
% p.Color(4) = 0.3 ;
% axis equal ;


X0 = [PT(:,1); PT(:,2); PT(:,3)] ;

% test = objfun(x1d, y1d, z1d, q, band, cpx_band, cpy_band, cpz_band, X0)

Xopt = fminsearch(@(X) objfun(x1d, y1d, z1d, q, band, cpx_band, cpy_band, cpz_band, X), X0) ;

numpt = length(Xopt)/3 ;
PT = [Xopt(1:numpt), Xopt(numpt+1:2*numpt), Xopt(2*numpt+1:3*numpt)] ;

figure(1) ; clf ;  hold on
patch('Faces', face, 'Vertices', vertex, 'FaceVertexCData', ...
    0.8*ones(size(face)), 'EdgeColor', 'none', 'FaceColor','flat') ;
plot3(PT(:,1), PT(:,2), PT(:,3), 'r-o', 'MarkerSize', 3, 'MarkerFaceColor', 'r')
axis equal ;


PT_new = PT ;

for k = 2:numpt-1
    Etemp = interp3_matrix(x1d, y1d, z1d, PT_new(k,1), PT_new(k,2), PT_new(k,3), q, band) ;
    PT_new(k,1) = Etemp*cpx_band ; PT_new(k,2) = Etemp*cpy_band ; PT_new(k,3) = Etemp*cpz_band ;
end

figure(1) ; clf ;  hold on
patch('Faces', face, 'Vertices', vertex, 'FaceVertexCData', ...
    0.8*ones(size(face)), 'EdgeColor', 'none', 'FaceColor','flat') ;
plot3(PT_new(:,1), PT_new(:,2), PT_new(:,3), 'k-o', 'MarkerSize', 3, 'MarkerFaceColor', 'k')
axis equal ;


%%

function objval = objfun(x1d, y1d, z1d, q, band, cpx_band, cpy_band, cpz_band, X)
% Get number of points
numpt = length(X)/3 ;
% Make input variable into (numpt)x3 matrix
PT = [X(1:numpt), X(numpt+1:2*numpt), X(2*numpt+1:3*numpt)] ;
% IBP for each point onto the surface (except the first and last, where are fixed)
for k = 2:numpt-1
    Etemp = interp3_matrix(x1d, y1d, z1d, PT(k,1), PT(k,2), PT(k,3), q, band) ;
    PT(k,1) = Etemp*cpx_band ; PT(k,2) = Etemp*cpy_band ; PT(k,3) = Etemp*cpz_band ;
end
% Compute discretized arclength
objval = 0 ;
for k = 2:numpt-1
objval = objval + (1/4*numpt)*norm(PT(k+1,:) - PT(k-1,:))^2 ;
end
end
