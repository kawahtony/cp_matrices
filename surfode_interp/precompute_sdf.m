% Compute closest points and signed distance function for triangular
% surfaces


clear ; clc; close all 

%% Load triangular mesh

% surface = 'trefoil_knot_tube' ;
surface = 'twisted_elliptic_torus' ;
% surface = 'wavy_torus' ;
% surface = 'double_torus_genus2' ;
% surface = 'genus3_torus' ;


[vertices, faces] = read_off(['data/', surface, '.off']) ;

figure(1) ; clf ;
trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3), 'FaceColor',[0.2 0.5 0.9], 'EdgeColor','none');
axis equal
camlight
lighting gouraud

temp = xlim ; xmin = temp(1) ; xmax = temp(2) ;
temp = ylim ; ymin = temp(1) ; ymax = temp(2) ;
temp = zlim ; zmin = temp(1) ; zmax = temp(2) ;


%% Cartesian mesh in the embedding space

dx = 0.02 ; % Grid size

dim = 3 ;  % Dimension of embedding space
q = 3 ; % Degree of interpolating polynomial

bw = 1 + sqrt(dim)*(q+1)/2 ; % Bandwidth coefficient

xmin = xmin - 0.5 ; xmax = xmax + 0.5 ;
ymin = ymin - 0.5 ; ymax = ymax + 0.5 ;
zmin = zmin - 0.5 ; zmax = zmax + 0.5 ;

x1d = (xmin:dx:xmax)' ; 
y1d = (ymin:dx:ymax)' ;
z1d = (zmin:dx:zmax)' ;

[xx, yy, zz] = meshgrid(x1d, y1d, z1d) ;


%% Compute sdf

[sdist, cpx, cpy, cpz, faceidx] = tri2sdf(xx, yy, zz, faces, vertices, bw*dx) ;


%% Banding

xg = xx(:) ;
yg = yy(:) ;
zg = zz(:) ;
cpxg = cpx(:) ;
cpyg = cpy(:) ;
cpzg = cpz(:) ;
sdistg = sdist(:) ;
faceidxg = faceidx(:) ;

band = find(abs(sdist) <= bw*dx) ;

% Store closest points in the band
xx_band = xg(band) ;
yy_band = yg(band) ;
zz_band = zg(band) ;
cpx_band = cpxg(band) ;
cpy_band = cpyg(band) ;
cpz_band = cpzg(band) ;
sdist_band = sdistg(band) ;
faceidx_band = faceidxg(band) ;



%% Compute curvature

% % The raw cotangent curvature can contain a mesh-scale alternating mode.
% % Ten local smoothing passes remove that mode while preserving the
% % large-scale curvature variation. Use 0 instead of 10 to inspect Hraw.
% [H, ~, ~, ~, Hraw] = mean_curvature_off(['data/', surface, '.off'], false, 10, 0.5) ;
% 
% 
% H_band = zeros(size(band)) ;
% for k = 1:length(band)
% 
%     faceidx_k = faceidx_band(k) ;
%     v1_idx = faces(faceidx_k, 1) ;
%     v2_idx = faces(faceidx_k, 2) ; 
%     v3_idx = faces(faceidx_k, 3) ;
% 
%     v1 = vertices(v1_idx, :)' ;
%     v2 = vertices(v2_idx, :)' ;
%     v3 = vertices(v3_idx, :)' ;
% 
%     cp_temp = [cpx_band(k); cpy_band(k); cpz_band(k)] ;
%     % Stable affine barycentric coordinates. This formulation explicitly
%     % enforces lambda(1) + lambda(2) + lambda(3) = 1.
%     edge12 = v2 - v1 ;
%     edge13 = v3 - v1 ;
%     rhs = cp_temp - v1 ;
%     gram = [dot(edge12, edge12), dot(edge12, edge13) ;
%             dot(edge12, edge13), dot(edge13, edge13)] ;
%     lambda23 = gram \ [dot(rhs, edge12); dot(rhs, edge13)] ;
%     lambda = [1 - sum(lambda23); lambda23] ;
%     H_band(k) = lambda(1)*H(v1_idx) + lambda(2)*H(v2_idx) + lambda(3)*H(v3_idx) ;
% end 


%% Linear operator matrices

Lmat = laplacian_3d_matrix(x1d, y1d, z1d, 2, band);
[Dxc, Dyc, Dzc] = firstderiv_cen2_3d_matrices(x1d, y1d, z1d, band) ;
Emat = interp3_matrix(x1d, y1d, z1d, cpx_band, cpy_band, cpz_band, q, band);


%% Heat smoothing of mean curvature

% H_band = Emat * H_band ; 
% 
% dt_smooth = 0.1 * dx^2 ; 
% N_smooth = 150 ;
% for k = 1:N_smooth
%     H_band = H_band + dt_smooth*(Lmat*H_band) ;
%     H_band = Emat * H_band ;
% end



%% Compute unit outward normal

nx = (xx_band - cpx_band).*sign(sdist_band) ; 
ny = (yy_band - cpy_band).*sign(sdist_band)  ;
nz = (zz_band - cpz_band).*sign(sdist_band)  ; 
nnorm = sqrt(nx.^2 + ny.^2 + nz.^2) ;

phi_x = Dxc*sdist_band ;
phi_y = Dyc*sdist_band ;
phi_z = Dzc*sdist_band ;
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


%% Save files

save_file = ['data/', surface, '_dx=', num2str(dx), '.mat'] ;
save(save_file)



%% Check sdf

% [faces_test, vertices_test] = isosurface(xx, yy, zz, reshape(sdist, size(xx)), -bw*dx) ;
% 
% figure(2) ; clf ; 
% % Plot surface
% patch('Faces', faces_test, 'Vertices', vertices_test, 'FaceVertexCData', ...
%     0.45*ones(size(faces_test)), 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceAlpha', 0.4) ;
% axis equal ; 


%% Safe space


% These surfaces are boring

% surface_file = 'data/harmonic_bumpy_sphere.off' ;
% surface_file = 'data/rounded_superellipsoid.off' ;