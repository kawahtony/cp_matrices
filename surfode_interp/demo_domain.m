
clear ; close all ;

azure = [0, 128, 255]/255 ;

%%

xvec = linspace(-1.5, 1.5, 50) ;
h = xvec(2) - xvec(1) ;

thvec = linspace(0, 2*pi, 100) ; thvec(end) = [] ;
Nth = length(thvec) ;

delta = h ; 
d = 2 ;
q = 3 ;
rad = delta + sqrt(d)*((q+1)/2)*h ;
rvec = linspace(1-delta, 1+delta, 5) ;

figure(1) ; clf ; hold on ;
axis equal ; axis off


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


R = 1 + 1.1*delta ; 
r = 1 - 1.1*delta ; 

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


%% Plot Omega_h

for I = 1:length(rvec)
    r = rvec(I) ;
    for J = 1:length(thvec)
        th = thvec(J) ;
        x = r*cos(th) ; y = r*sin(th) ;
        [i,j] = find_grid_box(xvec, x, y) ;
        tempset1 = [i-1, i, i+1, i+2] ;
        tempset2 = [j-1, j, j+1 ,j+2] ;
        for I2 = 1:4
            for J2 = 1:4
                plot(xvec(tempset1(I2)), xvec(tempset2(J2)), 'rx', 'MarkerSize', 7, 'LineWidth', 1.5)
            end
        end
    end
end

xlim([-0.3200, 1.6436])
ylim([0.1707, 1.2732])


%% Export the figure

exportgraphics(gcf, 'domain.pdf')


%% Helper functions


function [i,j] = find_grid_box(xvec, x, y)
    i = find(xvec < x, 1, 'last') ;
    j = find(xvec < y, 1, 'last') ;
end