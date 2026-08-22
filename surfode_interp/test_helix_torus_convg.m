azure = [0, 128, 255]/255 ; 


figure(1) ; clf ;

dtset = [0.02, 0.01, 0.005] ;


errset = zeros(size(dtset)) ;
for k = 1:length(dtset)
    dt = dtset(k) ;
    filename = ['results/helix_torus_RK4_dt=', num2str(dt), '.mat'] ;
    load(filename) ;
    err = norm(y_stored - y_exact_stored, 'inf') / norm(y_exact_stored, 'inf') ;
    errset(k) = err ;
end


loglog(dtset, errset, 'o', 'MarkerSize', 10, 'MarkerFaceColor', azure, 'Color', azure) ; hold on ;
loglog(dtset, dtset.^4, '-')