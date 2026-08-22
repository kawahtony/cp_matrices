
clear ; close all ;

azure = [0, 128, 255]/255 ;

%%

delta = 0.0612 ;

h = 0.05/2/sqrt(2) ;
d = 2 ;
q = 3 ;
rad = delta + sqrt(d)*((q+1)/2)*h ;

thvec = linspace(0, 2*pi, 100);
thvec(end) = [];   % avoid duplicate point at 2*pi
Nth = length(thvec) ;

figure(1) ; clf ;
hold on ;


%% Plot Omega_0


R = 1 + rad ; 
r = 1 - rad ; 

x_outer = R*cos(thvec);
y_outer = R*sin(thvec);

x_inner = r*cos(thvec);
y_inner = r*sin(thvec);


% Vertices: first outer circle, then inner circle
V = [x_outer(:), y_outer(:);
     x_inner(:), y_inner(:)];

% Faces: each quadrilateral connects two outer points and two inner points
F = zeros(Nth, 4);

for j = 1:Nth
    jp = mod(j, Nth) + 1;   % next index, periodic

    F(j,:) = [j, jp, Nth+jp, Nth+j];
end


patch('Vertices', V, ...
      'Faces', F, ...
      'FaceColor', azure, ...
      'EdgeColor', 'none', ...
      'FaceAlpha', 0.3);


%% Plot Omega


R = 1 + delta ; 
r = 1 - delta ; 

x_outer = R*cos(thvec);
y_outer = R*sin(thvec);

x_inner = r*cos(thvec);
y_inner = r*sin(thvec);


% Vertices: first outer circle, then inner circle
V = [x_outer(:), y_outer(:);
     x_inner(:), y_inner(:)];

% Faces: each quadrilateral connects two outer points and two inner points
F = zeros(Nth, 4);

for j = 1:Nth
    jp = mod(j, Nth) + 1;   % next index, periodic

    F(j,:) = [j, jp, Nth+jp, Nth+j];
end


patch('Vertices', V, ...
      'Faces', F, ...
      'FaceColor', azure, ...
      'EdgeColor', 'none', ...
      'FaceAlpha', 0.7);

%%


th = pi/4 ;

plot(exp(1i*linspace(0, 2*pi, 300)), 'k-', 'LineWidth', 3)

plot((1 + delta)*cos(th), (1 + delta)*sin(th), 'ko', 'MarkerSize', 12, 'LineWidth', 1.5)

xb = (1 + delta)*cos(th) - h ;
yb = (1 + delta)*sin(th) - h ;

for i = 0:3
    for j = 0:3
        plot(xb + i*h, yb + j*h, 'rx', 'LineWidth', 1.5, 'MarkerSize', 10)
    end
end

plot(xb + i*h, yb + j*h, 'ks', 'LineWidth', 1.5, 'MarkerSize', 12)


axis equal ; axis off ;
xlim([0.6474    0.8906])
ylim([0.6907    0.8444])


%% Export

exportgraphics(gcf, 'tube_radius.pdf')