close all ; clc ; clear ;  

azure = [0, 128, 255]/255 ; 
coral = [255, 127, 80]/255 ;
jade = [0, 168, 107]/255 ;

dxset = [0.08, 0.04, 0.02, 0.01] ;


%% Linear interpolation


figure(1) ; clf ;

errset = zeros(size(dxset)) ;
for k = 1:length(dxset)
    dx = dxset(k) ;
    filename = ['results/helix_torus/RK2_q=1_dx=', num2str(dx)] ;
    load(filename) ;
    diff = sqrt(  (x_stored(1,:) - x_exact_stored(1,:)).^2 + (x_stored(2,:) - x_exact_stored(2,:)).^2 ...
                + (x_stored(3,:) - x_exact_stored(3,:)).^2 ) ;   
    errset(k) = max(diff) ;
end

loglog(dxset, errset, 'o', 'MarkerSize', 10, 'MarkerFaceColor', azure, 'Color', azure, ...
       'DisplayName', 'RK2 with $q=1$') ; hold on ;

errset = zeros(size(dxset)) ;
for k = 1:length(dxset)
    dx = dxset(k) ;
    filename = ['results/helix_torus/RK3_q=1_dx=', num2str(dx)] ;
    load(filename) ;
    diff = sqrt(  (x_stored(1,:) - x_exact_stored(1,:)).^2 + (x_stored(2,:) - x_exact_stored(2,:)).^2 ...
                + (x_stored(3,:) - x_exact_stored(3,:)).^2 ) ;   
    errset(k) = max(diff) ;
end

loglog(dxset, errset, 's', 'MarkerSize', 10, 'MarkerFaceColor', coral, 'Color', coral, ...
       'DisplayName', 'RK3 with $q=1$') ; hold on ;

loglog(dxset, dxset.^2, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Reference: $h^2$')


xticks([0.01 0.02, 0.04, 0.08]) ; 


set(gca, 'FontSize', 20) 
set(gcf, 'defaultTextInterpreter', 'Latex')
set(gca, 'TickLabelInterpreter', 'Latex')
set(gcf, 'Position', [109   304   512   470])
leg = legend ; 
set(leg, 'Box', 'off', 'Interpreter', 'Latex', ...
    'Position', [0.1941    0.6063    0.4266    0.2720], 'FontSize', 14)
grid on ; 
xlim([0.008, 0.088])
ylim([7e-5, 2e-1])

xlabel('$h$') ; ylabel('Error')
% 
% exportgraphics(gcf, 'helix_torus_convg1.pdf', 'Resolution', 600)


%% Cubic interpolation


figure(2) ; clf ;

errset = zeros(size(dxset)) ;
for k = 1:length(dxset)
    dx = dxset(k) ;
    filename = ['results/helix_torus/RK2_q=3_dx=', num2str(dx)] ;
    load(filename) ;
    diff = sqrt(  (x_stored(1,:) - x_exact_stored(1,:)).^2 + (x_stored(2,:) - x_exact_stored(2,:)).^2 ...
                + (x_stored(3,:) - x_exact_stored(3,:)).^2 ) ;  
    errset(k) = max(diff) ;
end

loglog(dxset, errset, 'o', 'MarkerSize', 10, 'MarkerFaceColor', azure, 'Color', azure, ...
       'DisplayName', 'RK2 with $q=3$') ; hold on ;



errset = zeros(size(dxset)) ;
for k = 1:length(dxset)
    dx = dxset(k) ;
    filename = ['results/helix_torus/RK3_q=3_dx=', num2str(dx)] ;
    load(filename) ;
    diff = sqrt(  (x_stored(1,:) - x_exact_stored(1,:)).^2 + (x_stored(2,:) - x_exact_stored(2,:)).^2 ...
                + (x_stored(3,:) - x_exact_stored(3,:)).^2 ) ;  
    errset(k) = max(diff) ;
end

loglog(dxset, errset, 's', 'MarkerSize', 10, 'MarkerFaceColor', coral, 'Color', coral, ...
       'DisplayName', 'RK3 with $q=3$') ; hold on ;

errset = zeros(size(dxset)) ;
for k = 1:length(dxset)
    dx = dxset(k) ;
    filename = ['results/helix_torus/RK4_q=3_dx=', num2str(dx)] ;
    load(filename) ;
    diff = sqrt(  (x_stored(1,:) - x_exact_stored(1,:)).^2 + (x_stored(2,:) - x_exact_stored(2,:)).^2 ...
                + (x_stored(3,:) - x_exact_stored(3,:)).^2 ) ;   
    errset(k) = max(diff) ;
end

loglog(dxset, errset, 'd', 'MarkerSize', 10, 'MarkerFaceColor', jade, 'Color', jade, ...
       'DisplayName', 'RK4 with $q=3$') ; hold on ;

loglog(dxset, dxset.^2, '--', 'LineWidth', 1.5, 'Color', azure, 'DisplayName', 'Reference: $h^2$')

loglog(dxset, dxset.^3, '-.', 'LineWidth', 1.5, 'Color', coral, ...
       'DisplayName', 'Reference: $h^3$')

loglog(dxset, dxset.^4, ':', 'LineWidth', 1.5, 'Color', jade, ...
       'DisplayName', 'Reference: $h^4$')


xticks([0.01 0.02, 0.04, 0.08]) ; 


set(gca, 'FontSize', 20) 
set(gcf, 'defaultTextInterpreter', 'Latex')
set(gca, 'TickLabelInterpreter', 'Latex')
set(gcf, 'Position', [622   303   510   470])
leg = legend ; 
set(leg, 'Box', 'off', 'Interpreter', 'Latex', 'FontSize', 14, ...
    'NumColumns', 2, 'Position', [0.2732    0.1663    0.6260    0.1459]) 
grid on ; 
xlim([0.008, 0.088])
ylim([7e-11, 9e-2])
yticks([1e-10, 1e-8 1e-6 1e-4 1e-2])

xlabel('$h$') ; ylabel('Error')

% exportgraphics(gcf, 'helix_torus_convg2.pdf', 'Resolution', 600)
