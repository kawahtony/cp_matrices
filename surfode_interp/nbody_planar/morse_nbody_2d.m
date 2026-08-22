%% morse_nbody_2d.m
% N-particle 2D N-body simulation with a Morse interaction potential
% (short-range repulsion, long-range attraction), periodic boundary
% conditions, integrated with ode45, animated live. Same style/structure
% as coulomb_nbody_2d.m -- only the interaction law is different.
%
% Physics
%   Pairwise Morse potential (D'Orsogna-style swarm/aggregation potential):
%     U(r) = Cr*exp(-r/lr) - Ca*exp(-r/la)
%   Radial force magnitude (positive = repulsive, pushes particles apart):
%     F(r) = -dU/dr = (Cr/lr)*exp(-r/lr) - (Ca/la)*exp(-r/la)
%   Equation of motion:
%     m_i * d^2 r_i/dt^2 = sum_{j~=i} F(r_ij) * (r_i - r_j)/r_ij
%
%   With lr < la (repulsion shorter-range than attraction) and
%   Cr/lr > Ca/la (repulsion stronger at contact), F(r) > 0 (repulsive)
%   at small r and F(r) < 0 (attractive) at large r, i.e. exactly the
%   "short-range repulsion, long-range attraction" behavior requested.
%   There is a stable equilibrium separation r_eq where F(r_eq)=0:
%     r_eq = [ln(Ca/la) - ln(Cr/lr)] / (1/la - 1/lr)
%   Unlike the Coulomb 1/r^2 force, Morse force is finite at r=0, so no
%   softening is needed to avoid a singularity (only a tiny eps2 to guard
%   the r=0 divide when forming the unit vector).
%
% Periodic boundary conditions on a square box [0,L] x [0,L]:
%   - Pairwise separations use the minimum-image convention (nearest
%     periodic image), same approach as coulomb_nbody_2d.m.
%   - Positions are integrated unwrapped (smooth trajectories); wrapped
%     into the box only for plotting.
%
% Animation is driven by an ODE OutputFcn so the figure opens and starts
% updating the moment integration begins (see coulomb_nbody_2d.m notes on
% why: waiting until ode45 returns to build the figure can look like the
% script "did nothing" for a long time).

clear; clc; close all;

%% ---------------- Parameters ----------------
N   = 100;               % number of particles
L   = 20;               % periodic box size (domain is [0,L] x [0,L])

% Morse potential parameters: short-range repulsion, long-range attraction
Cr = 1.0;    lr = 0.5;   % repulsion strength & length scale (short range)
Ca = 0.5;    la = 2.0;   % attraction strength & length scale (long range)
eps2 = 1e-6;             % tiny softening, only to guard r=0 division

r_eq = (log(Ca/la) - log(Cr/lr)) / (1/la - 1/lr);
fprintf('Stable pairwise equilibrium separation r_eq = %.4f\n', r_eq);

rng(3);                                   % reproducible initial conditions
m = ones(N,1);                            % equal unit masses

P0 = L * rand(N,2);                       % random initial positions in box
V0 = 0.2 * (rand(N,2) - 0.5);             % small random initial velocities
y0 = [P0(:); V0(:)];                      % state = [x(:); y(:); vx(:); vy(:)]

tspan = 0:0.01:100;                        % output times for animation
opts  = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',0.1, ...
    'OutputFcn', @(t,y,flag) morseAnimOutputFcn(t,y,flag,N,L));

%% ---------------- Integrate (animates live via OutputFcn) ----------------
odefun = @(t,y) morsePeriodicODE(t,y,N,m,L,Cr,lr,Ca,la,eps2);
fprintf('Integrating... a figure window should appear now and animate live.\n');
[t, Y] = ode45(odefun, tspan, y0, opts);
fprintf('Integration finished (%d time points).\n', numel(t));

X  = Y(:,1:N);          % x_i(t), columns = particles
Yp = Y(:,N+1:2*N);      % y_i(t)

%% ---------------- Sanity check: total energy conservation ----------------
% Morse forces are smooth (no singularity), so this is a genuine
% conservative system and total energy should be well conserved.
V0m = reshape(y0(2*N+1:end), N, 2);
E0 = 0.5*sum(m.*sum(V0m.^2,2)) + morsePotentialEnergy(P0, L, Cr, lr, Ca, la, eps2);
Vend = reshape(Y(end,2*N+1:end), N, 2);
Pend = mod([X(end,:)', Yp(end,:)'], L);
Eend = 0.5*sum(m.*sum(Vend.^2,2)) + morsePotentialEnergy(Pend, L, Cr, lr, Ca, la, eps2);
fprintf('Total energy: start = %.6g, end = %.6g (rel. change %.3g%%)\n', ...
    E0, Eend, 100*abs(Eend-E0)/max(abs(E0),eps));

%% ================= Local functions =================
function status = morseAnimOutputFcn(t, y, flag, N, L)
% ODE OutputFcn that draws the live animation as ode45 integrates.
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
            title(hAx,'2D Morse N-body simulation (periodic boundary)  t = 0');
            xlabel(hAx,'x'); ylabel(hAx,'y');

            cmap = lines(N);
            hTrail = gobjects(N,1); hPts = gobjects(N,1);
            histX = cell(N,1); histY = cell(N,1);
            for i = 1:N
                col = cmap(i,:);
                hTrail(i) = plot(hAx, NaN, NaN, '-', 'Color', [col 0.4], 'LineWidth', 1);
                hPts(i) = plot(hAx, Xw0(i), Yw0(i), 'o', 'MarkerFaceColor', col, ...
                    'MarkerEdgeColor','k', 'MarkerSize', 8);
                histX{i} = Xw0(i); histY{i} = Yw0(i);
            end
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
            title(hAx, sprintf('2D Morse N-body simulation (periodic boundary)  t = %.2f / %.0f', t(end), tEnd));
            drawnow limitrate;

        case 'done'
            % nothing to clean up
    end
end

function dydt = morsePeriodicODE(~, y, N, m, L, Cr, lr, Ca, la, eps2)
% State layout: y = [x(:); y(:); vx(:); vy(:)], each block length N.
    P = reshape(y(1:2*N), N, 2);
    V = reshape(y(2*N+1:end), N, 2);
    A = zeros(N,2);
    for i = 1:N
        diff = P(i,:) - P;                    % N x 2
        diff = diff - L*round(diff/L);        % minimum-image convention
        r = sqrt(sum(diff.^2, 2) + eps2);     % pairwise distance
        r(i) = Inf;                           % no self-force
        Fmag = (Cr/lr).*exp(-r/lr) - (Ca/la).*exp(-r/la);  % >0 repulsive, <0 attractive
        rhat = diff ./ r;                     % N x 2 unit vectors
        A(i,:) = sum(Fmag .* rhat, 1) / m(i);
    end
    dydt = [V(:); A(:)];
end

function U = morsePotentialEnergy(P, L, Cr, lr, Ca, la, eps2)
% U = sum_{i<j} [Cr*exp(-r_ij/lr) - Ca*exp(-r_ij/la)]  (minimum-image r_ij)
    N = size(P,1);
    U = 0;
    for i = 1:N-1
        diff = P(i,:) - P(i+1:end,:);
        diff = diff - L*round(diff/L);
        r = sqrt(sum(diff.^2,2) + eps2);
        U = U + sum(Cr*exp(-r/lr) - Ca*exp(-r/la));
    end
end
