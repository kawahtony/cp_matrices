function [K, V, F, Kraw] = gaussian_curvature_off( ...
    filename, smoothingIterations, smoothingStrength, makePlot)
%GAUSSIAN_CURVATURE_OFF Compute Gaussian curvature of a triangular OFF mesh.
%
%   K = GAUSSIAN_CURVATURE_OFF('mesh.off') computes vertex Gaussian
%   curvature using angle defect divided by barycentric vertex area.
%
%   K = GAUSSIAN_CURVATURE_OFF('mesh.off', N, LAMBDA) applies N
%   positive-cotangent smoothing passes with blend LAMBDA. Recommended
%   values are N = 5 to 15 and LAMBDA = 0.5.
%
%   K = GAUSSIAN_CURVATURE_OFF(..., MAKEPLOT) controls plotting. MAKEPLOT
%   defaults to true.
%
%   [K,V,F,KRAW] also returns the mesh and unfiltered Gaussian curvature.
%   Gaussian curvature has units of inverse length squared.

    if nargin < 1
        error('gaussian_curvature_off:MissingInput', ...
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

    [~, K, ~, ~, V, F, raw] = surface_curvatures_off( ...
        filename, smoothingIterations, smoothingStrength);
    Kraw = raw.gaussian;

    if makePlot
        figure('Color', 'w');
        plot_vertex_field(F, V, K);
        title('Gaussian curvature');
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
