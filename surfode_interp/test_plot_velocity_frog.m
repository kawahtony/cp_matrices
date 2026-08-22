

%% Load sdf of frog surface directly

clear ; clc ; close all

jade = [0, 168, 107]/255 ;
% 
load('data/frog_dx=0.025.mat')
% 
[face, vertex] = isosurface(xx, yy, zz, sdist, 0) ;


%%

skip = 10 ;

x1d_coarse = x1d(1:skip:end) ;
y1d_coarse = y1d(1:skip:end) ;
z1d_coarse = z1d(1:skip:end) ;

[xx_coarse, yy_coarse, zz_coarse] = meshgrid(x1d_coarse, y1d_coarse, z1d_coarse) ;

[sdist_coarse, cpx_coarse, cpy_coarse, cpz_coarse, ~] = tri2sdf(xx_coarse, yy_coarse, zz_coarse, faces, vertices, bw*dx) ;

xx_coarse = xx_coarse(:) ;
yy_coarse = yy_coarse(:) ;
zz_coarse = zz_coarse(:) ;
cpx_coarse = cpx_coarse(:) ;
cpy_coarse = cpy_coarse(:) ;
cpz_coarse = cpz_coarse(:) ;
sdist_coarse = sdist_coarse(:) ;

band_coarse = find(abs(sdist_coarse) <= bw*dx) ;

xx_band_coarse = xx_coarse(band_coarse) ;
yy_band_coarse = yy_coarse(band_coarse) ;
zz_band_coarse = zz_coarse(band_coarse) ;
cpx_band_coarse = cpx_coarse(band_coarse) ;
cpy_band_coarse = cpy_coarse(band_coarse) ;
cpz_band_coarse = cpz_coarse(band_coarse) ;
sdist_band_coarse = sdist_coarse(band_coarse) ;


%% Compute unit outward normal


[Dx, Dy, Dz] = firstderiv_cen2_3d_matrices(x1d_coarse, y1d_coarse, z1d_coarse, band_coarse) ;

nx = (xx_band_coarse - cpx_band_coarse).*sign(sdist_band_coarse) ; 
ny = (yy_band_coarse - cpy_band_coarse).*sign(sdist_band_coarse)  ;
nz = (zz_band_coarse - cpz_band_coarse).*sign(sdist_band_coarse)  ; 
nnorm = sqrt(nx.^2 + ny.^2 + nz.^2) ;

phi_x = Dx*sdist_band_coarse ;
phi_y = Dy*sdist_band_coarse ;
phi_z = Dz*sdist_band_coarse ;
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


%% Impose ambient velocity field

% f1_coarse = cpy_band_coarse ;
% f2_coarse = - cpx_band_coarse ;
% f3_coarse = -0 * cpz_band_coarse ;


f1_coarse = 0 * cpx_band_coarse ;
f2_coarse = -(cpz_band_coarse - 1.5) ;
f3_coarse = cpy_band_coarse ;


%% Projection onto tangent space

for k = 1:length(band_coarse)
    vtemp = [f1_coarse(k); f2_coarse(k); f3_coarse(k)] ;
    ntemp = [nx(k); ny(k); nz(k)] ;
    vtemp = vtemp - dot(vtemp, ntemp)*ntemp ;
    f1_coarse(k) = vtemp(1) ;
    f2_coarse(k) = vtemp(2) ;
    f3_coarse(k) = vtemp(3) ;
end


%%



% Normalize the vectors to unit length so every arrow is clearly visible
% (the raw rotational field grows with radius, hiding the small arrows).
mag = sqrt(f1_coarse.^2 + f2_coarse.^2 + f3_coarse.^2) ;
mag(mag == 0) = 1 ;   % avoid division by zero
u = f1_coarse ./ mag ;
v = f2_coarse ./ mag ;
w = f3_coarse ./ mag ;


%% Plot


figure(1) ; clf ; 

% Plot surface (more transparent so the arrows show through)
patch('Faces', face, 'Vertices', vertex, 'FaceVertexCData', ...
    0.45*ones(size(face)), 'EdgeColor', 'none', 'FaceColor', jade, 'FaceAlpha', 0.8) ;
hold on
axis equal


% Fixed-length arrows, thick, red, with prominent heads for visibility
% arrowLen = 4 * skip * dx ;   % scale to the coarse grid spacing
arrowLen = skip * dx ;
q = quiver3(cpx_band_coarse, cpy_band_coarse, cpz_band_coarse, ...
    u, v, w, 0) ;            % 0 disables MATLAB auto-scaling
q.AutoScale   = 'off' ;
q.UData = arrowLen * u ;
q.VData = arrowLen * v ;
q.WData = arrowLen * w ;
q.Color       = [0.85 0.1 0.1] ;
q.LineWidth   = 0.5 ;
q.MaxHeadSize = 0.3 ;

axis equal
view(3) ;
camlight ; lighting gouraud ;
title('Velocity field on frog surface') ;


