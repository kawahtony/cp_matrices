function X = great_circle_points(npts, normal)
%GREAT_CIRCLE_POINTS Generate npts points on a great circle of the unit sphere.
%
%   X = great_circle_points(npts, normal)
%
%   Input:
%       npts   : number of points to generate
%       normal : 3-by-1 or 1-by-3 vector normal to the plane of the great circle
%
%   Output:
%       X      : npts-by-3 array, each row is a point on the great circle

    % Normalize the normal vector
    normal = normal(:);
    normal = normal / norm(normal);

    % Pick a vector not parallel to normal
    v = [1; 0; 0];
    if abs(dot(v, normal)) > 0.9
        v = [0; 1; 0];
    end

    % Construct an orthonormal basis {e1, e2} for the plane
    e1 = v - dot(v, normal) * normal;
    e1 = e1 / norm(e1);

    e2 = cross(normal, e1);
    e2 = e2 / norm(e2);

    % Parameter values
    theta = linspace(0, 2*pi, npts + 1);
    theta(end) = [];   % avoid duplicating the first point

    % Great circle parametrization
    X = zeros(npts, 3);
    for k = 1:npts
        X(k, :) = (cos(theta(k)) * e1 + sin(theta(k)) * e2).';
    end
end