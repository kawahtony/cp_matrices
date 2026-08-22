

numpt = 40 ;

PT = [linspace(-1, 1, numpt)', linspace(-1, 1, numpt)'] ;
% PT = [cos(linspace(0, pi/2, numpt))', sin(linspace(0, pi/2, numpt))'] ;

% PT(2:end-1,:) = PT(2:end-1,:) + 0.1*randn(size(PT(2:end-1,:))) ;

figure(1) ; clf ; 
plot(PT(:,1), PT(:,2), 'g-o', 'MarkerFaceColor', 'g', 'MarkerSize', 3) ; hold on ;

Xint0 = [PT(2:end-1,1); PT(2:end-1,2)] ;
options = optimset('MaxIter', 50000);
Xintopt = fminsearch(@(Xint) objfun(numpt,PT(1,:), PT(end,:), Xint), Xint0, options) ;

PT_opt = [Xintopt(1:numpt-2), Xintopt(numpt-1:2*numpt-4)] ;

figure(1) ; 
plot(PT_opt(:,1), PT_opt(:,2), 'k-o', 'MarkerFaceColor', 'k', 'MarkerSize', 3)
axis equal ;


%%
% 
function objval = objfun(numpt, p, q, Xint)
    % Make input variable into (numpt)x2 matrix
    PT = [Xint(1:numpt-2), Xint(numpt-1:2*numpt-4)] ;
    PT = [p; PT; q] ;    
    objval = 0 ;
    % Energy from arclength
    for k = 1:numpt-1
        objval = objval + (numpt/2)*(norm(PT(k+1,:) - PT(k,:)))^2 ;
    end
    % Energy from equal chord length
    mu = 5; 
    for k = 2:numpt-1
        objval = objval + mu*(norm(PT(k+1,:) - PT(k,:))^2 - norm(PT(k,:) - PT(k-1,:))^2 ) ;
    end
end
