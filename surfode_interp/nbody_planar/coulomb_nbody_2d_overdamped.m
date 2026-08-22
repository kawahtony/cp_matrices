%% coulomb_nbody_2d_overdamped.m
% N-particle 2D overdamped Coulomb (electrostatic) N-body simulation with
% periodic boundary conditions, integrated with ode45, and animated live.
%
% Physics (first-order / overdamped limit: inertia negligible vs drag)
%   gamma_i * dr_i/dt = F_i = sum_{j~=i} k * q_i * q_j * (r_i - r_j) / |r_i-r_j|^3
%   =>  dr_i/dt = F_i / gamma_i        ("Coulomb velocity")
%
% This is a gradient-descent flow on the pairwise Coulomb potential energy
% U = sum_{i<j} k*q_i*q_j/|r_i-r_j|, so U(t) is monotonically non-increasing
% (verified numerically below) -- unlike the inertial (2nd-order) version,
% there is no energy conservation to check here; monotonic U decrease is
% the correct sanity check instead.
%
% Periodic boundary conditions on a square box [0,L] x [0,L]:
%   - Pairwise separations use the minimum-image convention (each pair
%     interacts through whichever periodic image of particle j is closest
%     to particle i). This is the standard cheap approximation used in
%     teaching-level periodic N-body/MD codes; it is NOT a full Ewald
%     summation, so it is not exact for a truly infinite lattice sum, but
%     it is the correct leading-order periodic force.
%   - Positions are integrated unwrapped (so trajectories are smooth);
%     they are wrapped into the box only for plotting.
%
% A softening length eps is added in quadrature to the distance to avoid
% the force singularity when two particles get arbitrarily close.

clear; clc; close all;

%% ---------------- Parameters ----------------
N        = 20;              % number of particles
L        = 20;             % periodic box size (domain is [0,L] x [0,L])
k        = 1;              % Coulomb constant (normalized units)
eps_soft = 0.15;           % softening length
eps2     = eps_soft^2;
gamma    = 1;              % drag/friction coefficient (same for all particles)

rng(2);                                   % reproducible initial conditions
q = ones(N,1);                            % charges (all +1; edit for mixed signs)

P0 = L * rand(N,2);                       % random initial positions in box
y0 = P0(:);                               % state = [x(:); y(:)]  (no velocities)

tspan = 0:0.01:100;                        % output times for animation

% See coulomb_nbody_2d.m for why we animate live via an OutputFcn rather
% than after ode45 returns: with tight tolerances + a long tspan, nothing
% would appear on screen until the whole solve finished, which can look
% like the script "did nothing." Using OutputFcn opens the figure at the
% start of the solve and updates it as integration proceeds.
opts  = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',0.2, ...
    'OutputFcn', @(t,y,flag) coulombAnimOutputFcn(t,y,flag,N,L,q));

%% ---------------- Integrate (animates live via OutputFcn) ----------------
odefun = @(t,y) coulombOverdampedODE(t,y,N,k,q,gamma,L,eps2);
fprintf('Integrating... a figure window should appear now and animate live.\n');
[t, Y] = ode45(odefun, tspan, y0, opts);
fprintf('Integration finished (%d time points).\n', numel(t));

X  = Y(:,1:N);          % x_i(t), columns = particles
Yp = Y(:,N+1:2*N);      % y_i(t)

%% ---------------- Sanity check: potential energy should decrease ----------------
U0   = periodicPotentialEnergy(P0, q, k, L, eps2);
Pend = mod([X(end,:)', Yp(end,:)'], L);
Uend = periodicPotentialEnergy(Pend, q, k, L, eps2);
fprintf('Potential energy: start = %.6g, end = %.6g (should decrease for overdamped flow)\n', U0, Uend);

%% ================= Local functions =================
function status = coulombAnimOutputFcn(t, y, flag, N, L, q)
% ODE OutputFcn that draws the live animation as ode45 integrates, so a
% figure appears immediately (at 'init') instead of only after the whole
% solve finishes. Note: this state vector has NO velocity block (first
% order system), so y is just [x(:); y(:)], length 2N.
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
            title(hAx,'2D overdamped Coulomb N-body simulation (periodic boundary)  t = 0');
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
            title(hAx, sprintf('2D overdamped Coulomb N-body simulation (periodic boundary)  t = %.2f / %.0f', t(end), tEnd));
            drawnow limitrate;

        case 'done'
            % nothing to clean up
    end
end

function dydt = coulombOverdampedODE(~, y, N, k, q, gamma, L, eps2)
% First-order state: y = [x(:); y(:)], length 2N.
% dr_i/dt = F_i / gamma_i
    P = reshape(y, N, 2);
    dPdt = zeros(N,2);
    for i = 1:N
        diff = P(i,:) - P;                    % N x 2
        diff = diff - L*round(diff/L);        % minimum-image convention
        distSq = sum(diff.^2, 2) + eps2;      % softened squared distance
        distCube = distSq.^1.5;
        distCube(i) = Inf;                    % no self-force
        coeff = k * q(i) * q ./ distCube;     % N x 1
        F_i = sum(coeff .* diff, 1);          % 1 x 2 net Coulomb force
        dPdt(i,:) = F_i / gamma;              % overdamped: velocity = F/gamma
    end
    dydt = dPdt(:);
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
