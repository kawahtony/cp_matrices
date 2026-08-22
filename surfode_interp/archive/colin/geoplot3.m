%%
% Data from Tim
A = load("geodesic_data2.txt") ;

beta = 2 ;

figure(1); clf;
[x,y,z] = ellipsoid(0,0,0,1,1,beta, 64);
surf(x,y,z) %, 'facealpha', 0.5)
hold on;

lw = 'linewidth';
plot3(A(:,1), A(:,2), A(:,3), lw, 3, 'DisplayName', 'Tim')
plot3(A(end,1), A(end,2), A(end,3), 'bx', lw, 3)

% load('mytest.mat')
% plot3(x_stored(:,1), x_stored(:,2), x_stored(:,3), 'r-', lw, 2) 

axis equal
xlabel('x'); ylabel('y'); zlabel('z')

colormap gray


K = 200
Xr = [2.4*rand(K,1)-1.2  2.4*rand(K,1)-1.2  4.4*rand(K,1)-2.2];
[cpx,cpy,cpz,sd] = cpEllipsoid(Xr(:,1), Xr(:,2), Xr(:,3), [beta 1], [],'z');
for j = 1:K
  if sd(j) > 0
    plot3(cpx(j), cpy(j), cpz(j), 'co')
    plot3([Xr(j,1) cpx(j)],[Xr(j,2) cpy(j)], [Xr(j,3) cpz(j)], 'c-')
  end
end



x0 = A(1,:)
v0 = A(2,:) - A(1,:)  % not quite right
v0 = v0 ./ norm(v0)

x = x0;
v = v0;

%% problem
% x'' = 0
% or
% x' = v
% v' = 0

f = @(t,u) [u(4) u(5) u(6) 0 0 0]

%% another problem, delay to ccw planar
%kparam = 1
%f = @(t,u) [u(4) u(5) u(6) -kparam*([u(4) u(5) u(6)] - [-u(2) u(1) 0])]

X = [x0 v0 norm(v0)]
numsteps = 1000
Tf = 20;
dt = Tf / numsteps;
for k = 1:(2*numsteps)
  % vtar = [-x(2) x(1) 0];
  % v-vtar

  %xnew = x + dt*v;
  %vnew = v + dt*0;
  %vnew = v + dt*(-kparam*(v-vtar));

  unew = onestep_rk2(f, [x v], 0, dt);
  assert(length(unew) == 6)
  xnew = unew(1:3);
  vnew = unew(4:6);

  %x = xnew;
  %v = vnew;
  [cpx,cpy,cpz,sd] = cpEllipsoid(xnew(1),xnew(2),xnew(3), [2 1], [], 'z');
  n = xnew - [cpx cpy cpz]; disp(norm(n))

  %% this is interesting, and it works...
  % v = [cpx cpy cpz] - x;
  % v = v ./ norm(v)

  %% project the v onto the current tangent plane
  v = vnew - (dot(vnew,n) / dot(n,n))*n;
  % unit length matches the Tf from arclength param code
  v = v ./ norm(v);

  x = [cpx cpy cpz];
  X(k + 1, :) = [x v norm(v)];

  %% Half-way through, swap direction
  if k == numsteps
    xe = X(end, 1:3);
    % ve = -X(end, 4:6);
    v = -v;
  end
end


plot3(X(1:numsteps,1), X(1:numsteps,2), X(1:numsteps,3), 'r-', lw, 4)
% I = (numsteps+1):(2*numsteps);
% plot3(X(I,1), X(I,2), X(I,3), 'g:', lw, 4)
% plot3(xe(1), xe(2), xe(3), 'ro', lw, 4)
% plot3(X(1,1), X(1,2), X(1,3), 'rs', lw, 4)
% plot3(X(end,1), X(end,2), X(end,3), 'go', lw, 4)
%plot3(X(1,1), X(1,2), X(1,3), 'gs', lw, 4)


err = x0 - X(end,1:3)
normerr = norm(err)

% N    rk2       rk4
% 1000 0.42696   0.42696
% 2000 0.19959   0.19959
% 4000 0.095739
% 8000 0.046794


% 0.04679360110126132
% 0.04679360101230069
% 0.04679360110126132


%%
% N
% 200  1.164
% 400  0.54768
% 800  0.25418

legend('S', 'Tim code', 'Tim code end', 'odesurf')

% legend('S', 'Tim code', 'Tim code end', 'odesurf', 'return')
