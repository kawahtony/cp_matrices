function [k1, k2, H, K, V, F, raw] = principal_curvatures_off( ...
    filename, smoothingIterations, smoothingStrength, makePlot)
%PRINCIPAL_CURVATURES_OFF Compute principal curvatures of a triangular mesh.
%
%   [K1,K2] = PRINCIPAL_CURVATURES_OFF('mesh.off') computes ordered vertex
%   principal curvatures K1 >= K2 from signed mean curvature H and Gaussian
%   curvature K:
%
%       K1,2 = H +/- sqrt(max(H^2-K,0)).
%
%   [...] = PRINCIPAL_CURVATURES_OFF('mesh.off', N, LAMBDA) applies N
%   positive-cotangent smoothing passes with blend LAMBDA to the underlying
%   mean-curvature normal and Gaussian-curvature fields. Recommended values
%   are N = 5 to 15 and LAMBDA = 0.5.
%
%   [...] = PRINCIPAL_CURVATURES_OFF(..., MAKEPLOT) controls plotting.
%   MAKEPLOT defaults to true.
%
%   [K1,K2,H,K,V,F,RAW] also returns signed mean curvature, Gaussian
%   curvature, the mesh, and all unfiltered quantities in RAW.
%
%   K1 and K2 are positive on convex regions when the consistently ordered
%   face normals point outward. Their units are inverse length.

    if nargin < 1
        error('principal_curvatures_off:MissingInput', ...
              'Supply the path to an OFF file.');
    end
    if nargin < 2
        smoothingIterations = 0;
    end
    if nargin < 3
        smoothingStrength = 0.5;
    end
    if nargin < 4
        makePlot = true;
    end
    validateattributes(makePlot, {'logical', 'numeric'}, {'scalar'}, ...
                       mfilename, 'makePlot', 4);
    makePlot = logical(makePlot);

    [H, K, k1, k2, V, F, raw] = surface_curvatures_off( ...
        filename, smoothingIterations, smoothingStrength);

    if makePlot
        figure('Color', 'w');
        layout = tiledlayout(1, 2, 'TileSpacing', 'compact', ...
                             'Padding', 'compact');

        nexttile;
        plot_vertex_field(F, V, k1);
        title('Maximum principal curvature, k_1');

        nexttile;
        plot_vertex_field(F, V, k2);
        title('Minimum principal curvature, k_2');

        title(layout, 'Principal curvatures');
    end
end


function plot_vertex_field(F, V, values)
%PLOT_VERTEX_FIELD Plot linearly interpolated vertex data on a mesh.
    patch('Faces', F, 'Vertices', V, ...
          'FaceVertexCData', values, ...
          'FaceColor', 'interp', ...
          'EdgeColor', 'none');
    axis equal off;
    view(3);
    colormap(turbo);
    colorbar;
    camlight('headlight');
    lighting gouraud;
end
