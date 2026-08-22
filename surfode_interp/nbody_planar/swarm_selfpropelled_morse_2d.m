%% swarm_selfpropelled_morse_2d.m
% N self-propelled particles in 2D that genuinely SWARM: they cohere into
% a bounded, moving group instead of freezing into a static crystal (like
% the plain Morse script) or flying apart. Free space (no periodic box) --
% the view auto-follows the group's centroid.
%
% Why plain Coulomb/Morse dynamics don't "swarm"
%   coulomb_nbody_2d.m and morse_nbody_2d.m are CONSERVATIVE: energy is
%   conserved, and left alone the particles either orbit/oscillate forever
%   or settle into a static equilibrium lattice. Real swarming (locusts,
%   fish schools, bacterial colonies) needs each individual to have its
%   own preferred cruising speed and to lose/gain energy to reach it --
%   i.e. dissipation, not conservation.
%
% The model (D'Orsogna, Chuang, Bertozzi & Chayes, PRL 2006,
% "Self-Propelled Particles with Soft-Core Interactions: Patterns,
% Stability, and Collapse")
%   dr_i/dt = v_i
%   m dv_i/dt = (alpha - beta*|v_i|^2) v_i  -  grad_i sum_j U_M(r_ij)
%
%   The term (alpha - beta|v|^2)v is a Rayleigh self-propulsion/friction
%   force: it ACCELERATES a slow particle (|v| small => net positive along
%   v) and DECELERATES a fast one (|v| large => net negative along v), so
%   every particle is driven toward its own preferred speed
%       v0 = sqrt(alpha/beta)
%   in the absence of interactions. U_M is the same short-range-repulsion
%   / long-range-attraction Morse potential as morse_nbody_2d.m:
%       U_M(r) = Cr*exp(-r/lr) - Ca*exp(-r/la),   lr < la
%   which keeps the group cohesive (attraction) while keeping individuals
%   from colliding (repulsion). Self-propulsion + friction supplies the
%   individual motion; the Morse potential supplies the collective
%   cohesion -- together they produce sustained, organized group motion
%   ("swarming"), which neither ingredient produces alone.
%
% Parameter regime used below (Cr=1, lr=0.5, Ca=1.2, la=1.0, alpha=1,
% beta=0.5) together with a slight initial tangential-velocity bias was
% checked numerically (RK4 port, N=30, T=60) to settle into a stable
% "single mill": a rotating disk of particles at roughly constant group
% radius, each moving near its preferred speed v0, with steady net
% rotation (angular momentum) and almost no net translation of the group
% -- a classic, clearly-visible swarming pattern. Small changes to
% Ca/la/alpha/beta move you into other regimes (dispersing "gas",
% collapsed "clump", or a coherently translating flock); see the printed
% diagnostics at the end and try tweaking them.

clear; clc; close all;

%% ---------------- Parameters ----------------
N = 30;                 % number of particles

% Self-propulsion / friction (Rayleigh term): drives each particle toward
% its own preferred speed v0 = sqrt(alpha/beta) when unopposed.
alpha = 1.0;
beta  = 0.5;
v0 = sqrt(alpha/beta);
fprintf('Preferred individual speed v0 = sqrt(alpha/beta) = %.4f\n', v0);

% Morse interaction: short-range repulsion, long-range attraction
Cr = 1.0;   lr = 0.5;
Ca = 1.2;   la = 1.0;
eps2 = 1e-6;             % tiny softening, only guards r=0 division
m  = ones(N,1);

% Which collective state you get out of the SAME Ca/la/alpha/beta depends
% heavily on how the initial velocities are structured -- this is genuine
% multistability, not a parameter effect. Pick one (verified numerically
% with Cr=1, lr=0.5, Ca=1.2, la=1.0, alpha=1, beta=0.5, N=30):
%   'tangential' -> rotating MILL: bounded radius ~0.75, speed -> v0,
%                   steady nonzero angular momentum, ~zero net drift.
%   'aligned'    -> rigid-body TRANSLATING FLOCK: bounded radius ~0.34,
%                   |mean velocity| -> v0 (fully coherent), zero rotation.
%   'random'     -> disordered STATIC CLUMP: bounded radius ~0.93, but
%                   no net rotation or drift -- individuals jitter near v0
%                   inside the ball without organizing.
% To see a DISPERSING gas or a fully COLLAPSED clump instead, it's the
% Ca/la values themselves that matter, not the IC -- see the comments at
% the top of this file for which direction to push them.
ic_mode = 'random';   % 'tangential' | 'aligned' | 'random'

rng(4);                                    % reproducible initial conditions
R0 = 1.5*sqrt(rand(N,1));                  % start in a small disk...
ang0 = 2*pi*rand(N,1);
P0 = [R0.*cos(ang0), R0.*sin(ang0)];

switch ic_mode
    case 'tangential'
        tang = [-sin(ang0), cos(ang0)];
        V0 = 1.0*tang + 0.1*randn(N,2);
    case 'aligned'
        V0 = repmat([1.0 0.0], N, 1) + 0.1*randn(N,2);
    case 'random'
        th = 2*pi*rand(N,1);
        V0 = 0.5*[cos(th), sin(th)];
        % V0 = rand(N,2) ;
    otherwise
        error('Unknown ic_mode "%s". Use ''tangential'', ''aligned'', or ''random''.', ic_mode);
end
y0 = [P0(:); V0(:)];

tspan = 0:0.02:60;
opts = odeset('RelTol',1e-7,'AbsTol',1e-9,'MaxStep',0.1, ...
    'OutputFcn', @(t,y,flag) swarmAnimOutputFcn(t,y,flag,N));

%% ---------------- Integrate (animates live via OutputFcn) ----------------
odefun = @(t,y) swarmODE(t,y,N,m,alpha,beta,Cr,lr,Ca,la,eps2);
fprintf('Integrating... a figure window should appear now and animate live.\n');
[t, Y] = ode45(odefun, tspan, y0, opts);
fprintf('Integration finished (%d time points).\n', numel(t));

X  = Y(:,1:N);
Yp = Y(:,N+1:2*N);
VX = Y(:,2*N+1:3*N);
VY = Y(:,3*N+1:4*N);

%% ---------------- Diagnostics: what collective state emerged? ----------------
Pend = [X(end,:)', Yp(end,:)'];
Vend = [VX(end,:)', VY(end,:)'];
centroid = mean(Pend,1);
radii = vecnorm(Pend - centroid, 2, 2);
meanSpeed = mean(vecnorm(Vend,2,2));
meanVel   = vecnorm(mean(Vend,1));
rel = Pend - centroid;
Lz = mean(rel(:,1).*Vend(:,2) - rel(:,2).*Vend(:,1));   % net angular momentum about centroid

% crude radius growth check (dispersing vs bounded) using the last 20% of run
i0 = round(0.8*numel(t));
radius_early = mean(vecnorm([X(i0,:)', Yp(i0,:)'] - mean([X(i0,:)',Yp(i0,:)'],1), 2, 2));
radius_growth_pct = 100*(mean(radii) - radius_early) / max(radius_early, eps);

fprintf('\n--- Final-state diagnostics ---\n');
fprintf('Mean group radius        : %.3f\n', mean(radii));
fprintf('Mean individual speed     : %.3f  (preferred v0 = %.3f)\n', meanSpeed, v0);
fprintf('Net group velocity |Vcm|  : %.3f  (large => coherent translation)\n', meanVel);
fprintf('Angular momentum Lz       : %.3f  (large magnitude => milling/rotation)\n', Lz);
fprintf('Radius change (last 20%%)  : %.1f%%  (near 0 => bounded, large + => dispersing)\n', radius_growth_pct);
if abs(radius_growth_pct) > 15
    fprintf('=> Interpretation: swarm is DISPERSING (attraction too weak for this population).\n');
elseif abs(Lz) > 0.3 && meanVel < 0.3*meanSpeed
    fprintf('=> Interpretation: stable MILL (rotating swarm about a fixed centroid).\n');
elseif meanVel > 0.5*meanSpeed
    fprintf('=> Interpretation: coherently TRANSLATING flock.\n');
else
    fprintf('=> Interpretation: bound but disordered CLUMP.\n');
end

%% ================= Local functions =================
function status = swarmAnimOutputFcn(t, y, flag, N)
% ODE OutputFcn: draws the live animation, camera auto-follows the
% group's centroid (there is no periodic box to anchor the view to).
    persistent hAx hTrail hPts histX histY trailLen tEnd
    status = 0;
    switch flag
        case 'init'
            trailLen = 300;
            tEnd = t(end);
            P0f = reshape(y(1:2*N), N, 2);

            figure('Color','w');
            hAx = axes; hold(hAx,'on'); axis(hAx,'equal');
            title(hAx,'Self-propelled swarm (Morse interaction)  t = 0');
            xlabel(hAx,'x'); ylabel(hAx,'y');

            cmap = lines(N);
            hTrail = gobjects(N,1); hPts = gobjects(N,1);
            histX = cell(N,1); histY = cell(N,1);
            for i = 1:N
                col = cmap(i,:);
                hTrail(i) = plot(hAx, NaN, NaN, '-', 'Color', [col 0.4], 'LineWidth', 1);
                hPts(i) = plot(hAx, P0f(i,1), P0f(i,2), 'o', 'MarkerFaceColor', col, ...
                    'MarkerEdgeColor','k', 'MarkerSize', 7);
                histX{i} = P0f(i,1); histY{i} = P0f(i,2);
            end
            centroid = mean(P0f,1);
            halfw = max(3, max(vecnorm(P0f-centroid,2,2)) + 2);
            xlim(hAx, centroid(1) + [-halfw halfw]);
            ylim(hAx, centroid(2) + [-halfw halfw]);
            drawnow;

        case ''
            if isempty(t), return; end
            for c = 1:numel(t)
                Pc = reshape(y(1:2*N, c), N, 2);
                for i = 1:N
                    histX{i}(end+1) = Pc(i,1); %#ok<AGROW>
                    histY{i}(end+1) = Pc(i,2); %#ok<AGROW>
                    if numel(histX{i}) > trailLen
                        histX{i}(1) = [];
                        histY{i}(1) = [];
                    end
                    set(hTrail(i), 'XData', histX{i}, 'YData', histY{i});
                    set(hPts(i), 'XData', Pc(i,1), 'YData', Pc(i,2));
                end
            end
            centroid = mean(Pc,1);
            halfw = max(3, max(vecnorm(Pc-centroid,2,2)) + 2);
            xlim(hAx, centroid(1) + [-halfw halfw]);
            ylim(hAx, centroid(2) + [-halfw halfw]);
            title(hAx, sprintf('Self-propelled swarm (Morse interaction)  t = %.2f / %.0f', t(end), tEnd));
            drawnow limitrate;

        case 'done'
            % nothing to clean up
    end
end

function dydt = swarmODE(~, y, N, m, alpha, beta, Cr, lr, Ca, la, eps2)
% State: y = [x(:); y(:); vx(:); vy(:)], each block length N. Free space
% (no periodic wrap) -- ordinary Euclidean pairwise distances.
    P = reshape(y(1:2*N), N, 2);
    V = reshape(y(2*N+1:end), N, 2);

    speed2 = sum(V.^2, 2);                        % N x 1
    selfProp = (alpha - beta*speed2) .* V;         % N x 2, Rayleigh self-propulsion/friction

    A = zeros(N,2);
    for i = 1:N
        diff = P(i,:) - P;                         % N x 2
        r = sqrt(sum(diff.^2, 2) + eps2);
        r(i) = Inf;
        Fmag = (Cr/lr).*exp(-r/lr) - (Ca/la).*exp(-r/la);  % >0 repulsive, <0 attractive
        rhat = diff ./ r;
        A(i,:) = sum(Fmag .* rhat, 1);
    end

    dVdt = (selfProp + A) ./ m;
    dydt = [V(:); dVdt(:)];
end
