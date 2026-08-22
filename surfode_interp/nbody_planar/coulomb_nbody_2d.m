%% coulomb_nbody_2d.m
% N-particle 2D Coulomb (electrostatic) N-body simulation with periodic
% boundary conditions, integrated with ode45, and animated live.
%
% Physics
%   m_i * d^2 r_i/dt^2 = sum_{j~=i} k * q_i * q_j * (r_i - r_j) / |r_i-r_j|^3
%
% Periodic boundary conditions on a square box [0,L] x [0,L]:
%   - Pairwise separations use the minimum-image convention (each pair
%     interacts through whichever periodic image of particle j is closest
%     to particle i). This is the standard cheap approximation used in
%     teaching-level periodic N-body/MD codes; it is NOT a full Ewald
%     summation, so it is not exact for a truly infinite lattice sum, but
%     it is the correct leading-order periodic force and is what people
%     mean by "Coulomb with periodic BCs" in this context.
%   - Positions are integrated unwrapped (so trajectories are smooth);
%     they are wrapped into the box only for plotting.
%
% A softening length eps is added in quadrature to the distance to avoid
% the force singularity when two particles get arbitrarily close.

clear; clc; close all;

%% ---------------- Parameters ----------------
N       = 20 ;            % number of particles
L       = 20;            % periodic box size (domain is [0,L] x [0,L])
k       = 1;              % Coulomb constant (normalized units)
eps_soft = 0.15;           % softening length
eps2    = eps_soft^2;

rng(2);                                   % reproducible initial conditions
% q = 2*(randi([0,1],N,1)) - 1;             % random charges, +-1
q = ones(N,1) ;
m = ones(N,1);                            % equal unit masses

P0 = L * rand(N,2);                       % random initial positions in box
V0 = 0.6 * (rand(N,2) - 0.5);             % small random initial velocities
y0 = [P0(:); V0(:)];                      % state = [x(:); y(:); vx(:); vy(:)]

tspan = 0:0.01:100;                        % output times for animation

% NOTE ON SPEED: RelTol/AbsTol this tight (1e-9/1e-11) combined with a
% long tspan and a chaotic repulsive many-body force can make ode45 take
% a very long time to finish BEFORE anything is drawn, if the figure is
% only built after ode45 returns. To avoid "nothing happens for minutes",
% this script instead animates live via an OutputFcn, so a window opens
% and starts updating immediately as ode45 integrates -- and we use more
% modest (still accurate) tolerances so it doesn't crawl.
opts  = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',0.1, ...
    'OutputFcn', @(t,y,flag) coulombAnimOutputFcn(t,y,flag,N,L,q));

%% ---------------- Integrate (animates live via OutputFcn) ----------------
odefun = @(t,y) coulombPeriodicODE(t,y,N,k,q,m,L,eps2);
fprintf('Integrating... a figure window should appear now and animate live.\n');
[t, Y] = ode45(odefun, tspan, y0, opts);
fprintf('Integration finished (%d time points).\n', numel(t));

X  = Y(:,1:N);          % x_i(t), columns = particles
Yp = Y(:,N+1:2*N);      % y_i(t)

%% ---------------- Sanity check: pairwise potential energy ----------------
% (Not conserved exactly since we didn't symplectically integrate, but
% ode45 with tight tolerances should keep drift small over the run.)
U0 = periodicPotentialEnergy(P0, q, k, L, eps2);
Pend = mod([X(end,:)', Yp(end,:)'], L);
Uend = periodicPotentialEnergy(Pend, q, k, L, eps2);
fprintf('Potential energy: start = %.6g, end = %.6g (rel. change %.3g%%)\n', ...
    U0, Uend, 100*abs(Uend-U0)/max(abs(U0),eps));

%% ================= Local functions =================
function status = coulombAnimOutputFcn(t, y, flag, N, L, q)
% ODE OutputFcn that draws the live animation as ode45 integrates, so a
% figure appears immediately (at 'init') instead of only after the whole
% solve finishes.
    persistent hAx hTrail hPts histX histY trailLen tEnd
    status = 0;
    switch flag
        case 'init'
            trailLen = 300;
            tEnd = t(end);
            P0f = reshape(y(1:2*N), N, 2);
            Xw0 = mod(P0f(:,1), L); Yw0 = mod(P0f(:,2), L);

            figure('Color','w');
            hAx = axes; hold(hAx,'on'); axis(hAx,'equal');
            xlim(hAx,[0 L]); ylim(hAx,[0 L]);
            rectangle('Position',[0 0 L L],'EdgeColor',[0.3 0.3 0.3],'LineWidth',1.2);
            title(hAx,'2D Coulomb N-body simulation (periodic boundary)  t = 0');
            xlabel(hAx,'x'); ylabel(hAx,'y');

            cmap = lines(N);
            hTrail = gobjects(N,1); hPts = gobjects(N,1);
            histX = cell(N,1); histY = cell(N,1);
            for i = 1:N
                col = cmap(i,:);
                hTrail(i) = plot(hAx, NaN, NaN, '-', 'Color', [col 0.4], 'LineWidth', 1);
                mk = 'o'; if q(i) < 0, mk = 's'; end
                hPts(i) = plot(hAx, Xw0(i), Yw0(i), mk, 'MarkerFaceColor', col, ...
                    'MarkerEdgeColor','k', 'MarkerSize', 8);
                histX{i} = Xw0(i); histY{i} = Yw0(i);
            end
            legend(hPts(1:min(N,2)), {'+ charge','- charge'}, 'Location','bestoutside');
            drawnow;

        case ''
            if isempty(t), return; end
            for c = 1:numel(t)
                Pc = reshape(y(1:2*N, c), N, 2);
                Xw = mod(Pc(:,1), L); Yw = mod(Pc(:,2), L);
                for i = 1:N
                    if abs(Xw(i)-histX{i}(end)) > L/2 || abs(Yw(i)-histY{i}(end)) > L/2
                        histX{i}(end+1) = NaN; %#ok<AGROW>
                        histY{i}(end+1) = NaN; %#ok<AGROW>
                    end
                    histX{i}(end+1) = Xw(i); %#ok<AGROW>
                    histY{i}(end+1) = Yw(i); %#ok<AGROW>
                    if numel(histX{i}) > trailLen
                        histX{i}(1:end-trailLen) = [];
                        histY{i}(1:end-trailLen) = [];
                    end
                    set(hTrail(i), 'XData', histX{i}, 'YData', histY{i});
                    set(hPts(i), 'XData', Xw(i), 'YData', Yw(i));
                end
            end
            title(hAx, sprintf('2D Coulomb N-body simulation (periodic boundary)  t = %.2f / %.0f', t(end), tEnd));
            drawnow limitrate;

        case 'done'
            % nothing to clean up
    end
end

function dydt = coulombPeriodicODE(~, y, N, k, q, m, L, eps2)
% State layout: y = [x(:); y(:); vx(:); vy(:)], each block length N.
    P = reshape(y(1:2*N), N, 2);
    V = reshape(y(2*N+1:end), N, 2);
    A = zeros(N,2);
    for i = 1:N
        diff = P(i,:) - P;                    % N x 2
        diff = diff - L*round(diff/L);        % minimum-image convention
        distSq = sum(diff.^2, 2) + eps2;      % softened squared distance
        distCube = distSq.^1.5;
        distCube(i) = Inf;                    % no self-force
        coeff = k * q(i) * q ./ distCube;     % N x 1
        A(i,:) = sum(coeff .* diff, 1) / m(i);
    end
    dydt = [V(:); A(:)];
end

function U = periodicPotentialEnergy(P, q, k, L, eps2)
% U = sum_{i<j} k*q_i*q_j / sqrt(|r_i-r_j|^2_min_image + eps2)
    N = size(P,1);
    U = 0;
    for i = 1:N-1
        diff = P(i,:) - P(i+1:end,:);
        diff = diff - L*round(diff/L);
        dist = sqrt(sum(diff.^2,2) + eps2);
        U = U + sum(k * q(i) * q(i+1:end) ./ dist);
    end
end
