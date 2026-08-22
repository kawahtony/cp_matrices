% 
% 
% CASE = 'descent' ;
% 
% load(['critical_temp_', CASE, '.mat']) ;
% 
% critical_stored_refined = zeros(size(critical_stored)) ;
% 
% for k = 1:size(critical_stored, 2)
% 
%     temp = critical_stored(:,k) ;
%     xc0 = temp(1) ; yc0 = temp(2) ; zc0 = temp(3) ;
% 
%     [~, ~, ~, ~, u0, v0] = cpWavyTorus(xc0, yc0, zc0) ;
% 
% 
%     temp = fsolve(@wavyTorusMeanCurvature, [u0; v0]) ;
%     u = temp(1) ; v = temp(2) ;
% 
% 
%     R  = 2.05 + 0.24*cos(5*u) ;
%     zc = 0.18*sin(5*u) ;
%     r = 0.52*(1 + 0.13*cos(3*v + 2*u)) ;
%     xc = (R + r*cos(v))*cos(u) ;
%     yc = (R + r*cos(v))*sin(u) ;
%     zc = zc + r*sin(v) ;
%     critical_stored_refined(:,k) = [xc; yc; zc] ;
% end
% 
% save(['critical_temp_', CASE, '_refined.mat'], 'critical_stored_refined')


% %% Load surface data
% 
% load('data/wavy_torus_dx=0.025.mat') ;
% 
% Lmat = laplacian_3d_matrix(x1d, y1d, z1d, 2, band);
% Emat = interp3_matrix(x1d, y1d, z1d, cpx_band, cpy_band, cpz_band, q, band);
% [Dxc, Dyc, Dzc] = firstderiv_cen2_3d_matrices(x1d, y1d, z1d, band) ;
% 
% Lcpx = Emat * (Lmat * cpx_band) ; 
% Lcpy = Emat * (Lmat * cpy_band) ; 
% Lcpz = Emat * (Lmat * cpz_band) ; 
% 
% H = zeros(size(band)) ;
% 
% for k = 1:length(band)
%     H(k) = -dot([Lcpx(k); Lcpy(k); Lcpz(k)], [nx(k); ny(k); nz(k)]) ;
% end
% 
% H = Emat * H ; 
% 
% 
% %% Initialize plot
% 
% xp = vertices(:,1) ;
% yp = vertices(:,2) ;
% zp = vertices(:,3) ;
% 
% Eplot = interp3_matrix(x1d, y1d, z1d, xp, yp, zp, q, band) ;
% Hplot = Eplot*H ;
% 
% 
% scale = 0.985 ;
% 
% figure(1) ; clf ; hold on ;
% trisurf(faces, scale*xp, scale*yp, scale*zp, Hplot, 'FaceAlpha', 1, 'EdgeColor', 'none') ;
% plot3(xc0, yc0, zc0, 'o', 'MarkerFaceColor', 'r', 'MarkerSize', 10) ;
% plot3(xc, yc, zc, 'o', 'MarkerFaceColor', 'k', 'MarkerSize', 10) ;
% 
% shading interp
% axis equal off  
% cb = colorbar ;
% set(cb, 'FontSize', 20)
% cb.Label.Interpreter = 'Latex' ;
% colormap cool
% view(0.2000, 58.3811) 
% 
% 


%%



CASE = 'ascent' ;

load(['critical_temp_', CASE, '.mat']) ;

load(['critical_temp_', CASE, '_refined.mat']) ;

rel_res = zeros(size(res_stored)) ;

for k = 1:length(critical_stored)    
    rel_res(k) = norm(critical_stored(:,k) - critical_stored_refined(:,k)) / norm(critical_stored_refined(:,k)) ;
end