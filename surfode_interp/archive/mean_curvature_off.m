function [H, V, F, Hn, Hraw] = mean_curvature_off( ...
    filename, useSignedCurvature, smoothingIterations, smoothingStrength)
%MEAN_CURVATURE_OFF Compute and plot mean curvature of a triangular OFF mesh.
%
%   H = MEAN_CURVATURE_OFF('mesh.off') reads a triangular OFF file,
%   computes unsigned mean curvature at every vertex using the cotangent
%   Laplace--Beltrami operator.
%
%   H = MEAN_CURVATURE_OFF('mesh.off', true) returns signed curvature.
%   Signed curvature is positive on convex regions when the face ordering
%   is consistent and the vertex normals point outward.
%
%   H = MEAN_CURVATURE_OFF('mesh.off', true, N, LAMBDA) applies N explicit
%   smoothing iterations to the mean-curvature normal field. Each iteration
%   blends a fraction LAMBDA of the current value with a positive-cotangent
%   weighted neighbor average. Recommended values are N = 5 to 15 and
%   LAMBDA = 0.5. Set N = 0 to return the unfiltered estimate.
%
%   [H,V,F,Hn,Hraw] also returns the vertices, triangular faces, filtered
%   mean-curvature normal vectors, and the unfiltered scalar curvature.
%   OFF indices are assumed to be zero-based.
%
%   Notes:
%   * The input must be an ASCII OFF file containing triangular faces.
%   * Curvature has units of inverse length.
%   * Smoothing reduces mesh-scale oscillation but does not make a
%     piecewise-linear interpolant C1 across triangle edges.
%   * Values at open boundaries should be interpreted cautiously.

    if nargin < 1
        error('mean_curvature_off:MissingInput', ...
              'Supply the path to an OFF file.');
    end
    if nargin < 2
        useSignedCurvature = false;
    end
    if nargin < 3
        smoothingIterations = 0;
    end
    if nargin < 4
        smoothingStrength = 0.5;
    end
    validateattributes(useSignedCurvature, {'logical', 'numeric'}, ...
                       {'scalar'}, mfilename, 'useSignedCurvature', 2);
    validateattributes(smoothingIterations, {'numeric'}, ...
                       {'scalar', 'real', 'finite', 'integer', 'nonnegative'}, ...
                       mfilename, 'smoothingIterations', 3);
    validateattributes(smoothingStrength, {'numeric'}, ...
                       {'scalar', 'real', 'finite', '>=', 0, '<=', 1}, ...
                       mfilename, 'smoothingStrength', 4);
    useSignedCurvature = logical(useSignedCurvature);

    [V, F] = read_triangular_off(filename);
    numberOfVertices = size(V, 1);

    % Triangle vertex positions and twice the triangle areas.
    p1 = V(F(:, 1), :);
    p2 = V(F(:, 2), :);
    p3 = V(F(:, 3), :);
    faceNormal = cross(p2 - p1, p3 - p1, 2);
    doubleArea = sqrt(sum(faceNormal.^2, 2));

    meshScale = norm(max(V, [], 1) - min(V, [], 1));
    areaTolerance = max(realmin, 1e-14 * meshScale^2);
    nondegenerate = doubleArea > areaTolerance;
    if any(~nondegenerate)
        warning('mean_curvature_off:DegenerateFaces', ...
                '%d degenerate face(s) were ignored in cotangent weights.', ...
                sum(~nondegenerate));
    end

    % Cotangent of the angle at each corner. A corner's cotangent is the
    % weight of the edge opposite that corner.
    cot1 = zeros(size(doubleArea));
    cot2 = zeros(size(doubleArea));
    cot3 = zeros(size(doubleArea));
    cot1(nondegenerate) = row_dot(p2(nondegenerate,:) - p1(nondegenerate,:), ...
                                  p3(nondegenerate,:) - p1(nondegenerate,:)) ...
                                  ./ doubleArea(nondegenerate);
    cot2(nondegenerate) = row_dot(p1(nondegenerate,:) - p2(nondegenerate,:), ...
                                  p3(nondegenerate,:) - p2(nondegenerate,:)) ...
                                  ./ doubleArea(nondegenerate);
    cot3(nondegenerate) = row_dot(p1(nondegenerate,:) - p3(nondegenerate,:), ...
                                  p2(nondegenerate,:) - p3(nondegenerate,:)) ...
                                  ./ doubleArea(nondegenerate);

    % Symmetric cotangent-weight matrix. Contributions from the two faces
    % incident on an interior edge are added automatically by sparse().
    i = [F(:,2); F(:,3); F(:,3); F(:,1); F(:,1); F(:,2)];
    j = [F(:,3); F(:,2); F(:,1); F(:,3); F(:,2); F(:,1)];
    s = [cot1; cot1; cot2; cot2; cot3; cot3];
    W = sparse(i, j, s, numberOfVertices, numberOfVertices);
    L = spdiags(sum(W, 2), 0, numberOfVertices, numberOfVertices) - W;

    % Barycentric area associated with each vertex.
    faceArea = 0.5 * doubleArea;
    vertexArea = accumarray(F(:), repmat(faceArea / 3, 3, 1), ...
                            [numberOfVertices, 1], @sum, 0);

    % L*V/(4*A) approximates the mean-curvature normal H*n. The factor 4
    % combines the cotangent Laplacian's 1/(2*A) normalization with
    % Delta_surface(x) = -2*H*n (up to the normal/sign convention).
    Hn = nan(numberOfVertices, 3);
    validVertex = vertexArea > areaTolerance;
    Hn(validVertex, :) = bsxfun(@rdivide, L(validVertex,:) * V, ...
                                4 * vertexArea(validVertex));

    % Area-weighted vertex normals, inherited from face orientation.
    vertexNormal = zeros(numberOfVertices, 3);
    repeatedNormals = [faceNormal; faceNormal; faceNormal];
    for coordinate = 1:3
        vertexNormal(:, coordinate) = accumarray( ...
            F(:), repeatedNormals(:, coordinate), ...
            [numberOfVertices, 1], @sum, 0);
    end
    normalLength = sqrt(sum(vertexNormal.^2, 2));
    hasNormal = normalLength > 0;
    vertexNormal(hasNormal,:) = bsxfun(@rdivide, ...
        vertexNormal(hasNormal,:), normalLength(hasNormal));

    if useSignedCurvature
        Hraw = sum(Hn .* vertexNormal, 2);
    else
        Hraw = sqrt(sum(Hn.^2, 2));
    end
    Hraw(~hasNormal) = NaN;

    % Filter the vector field rather than only the scalar magnitude. This
    % keeps Hn and H consistent and damps tangential discretization noise.
    if smoothingIterations > 0 && smoothingStrength > 0
        Hn = smooth_cotangent_field( ...
            Hn, W, smoothingIterations, smoothingStrength);
    end

    if useSignedCurvature
        H = sum(Hn .* vertexNormal, 2);
        plotTitle = 'Signed mean curvature';
    else
        H = sqrt(sum(Hn.^2, 2));
        plotTitle = 'Mean curvature magnitude';
    end
    H(~hasNormal) = NaN;

    % Plot vertex values with linear interpolation over each triangle.
    % figure(1); clf ;
    % patch('Faces', F, 'Vertices', V, ...
    %       'FaceVertexCData', H, ...
    %       'FaceColor', 'interp', ...
    %       'EdgeColor', 'none');
    % axis equal 
    % axis off;
    % view(3);
    % colormap(turbo);
    % colorbar;
    % title(plotTitle);
    % camlight('headlight');
    % lighting gouraud;
end


function field = smooth_cotangent_field(field, W, iterations, lambda)
%SMOOTH_COTANGENT_FIELD Stable local averaging of a vertex vector field.
%   Negative cotangent weights are discarded only for smoothing. The
%   curvature operator above still uses the complete cotangent matrix.

    positiveW = spfun(@(weight) max(weight, 0), W);
    finiteVertex = all(isfinite(field), 2);

    for iteration = 1:iterations %#ok<NASGU>
        safeField = field;
        safeField(~finiteVertex,:) = 0;

        availableWeight = positiveW * double(finiteVertex);
        canSmooth = finiteVertex & availableWeight > 0;
        neighborSum = positiveW * safeField;
        neighborAverage = field;
        neighborAverage(canSmooth,:) = bsxfun( ...
            @rdivide, neighborSum(canSmooth,:), availableWeight(canSmooth));

        field(canSmooth,:) = (1 - lambda) * field(canSmooth,:) + ...
                             lambda * neighborAverage(canSmooth,:);
    end
end


function d = row_dot(a, b)
%ROW_DOT Dot product of corresponding matrix rows.
    d = sum(a .* b, 2);
end


function [V, F] = read_triangular_off(filename)
%READ_TRIANGULAR_OFF Read vertices and triangular faces from an ASCII OFF file.
    fileId = fopen(filename, 'r');
    if fileId == -1
        error('mean_curvature_off:CannotOpenFile', ...
              'Could not open OFF file: %s', filename);
    end
    closeFile = onCleanup(@() fclose(fileId)); %#ok<NASGU>

    contents = fread(fileId, '*char')';
    % Remove comments beginning with #, then tokenize all remaining text.
    contents = regexprep(contents, '#[^\r\n]*', ' ');
    tokens = regexp(contents, '\S+', 'match');

    if isempty(tokens) || ~strcmpi(tokens{1}, 'OFF')
        error('mean_curvature_off:InvalidHeader', ...
              'The file must begin with the ASCII OFF header.');
    end
    if numel(tokens) < 4
        error('mean_curvature_off:IncompleteFile', ...
              'The OFF file does not contain valid mesh counts.');
    end

    numberOfVertices = str2double(tokens{2});
    numberOfFaces = str2double(tokens{3});
    if ~isfinite(numberOfVertices) || ~isfinite(numberOfFaces) || ...
       numberOfVertices < 1 || numberOfFaces < 1 || ...
       numberOfVertices ~= floor(numberOfVertices) || ...
       numberOfFaces ~= floor(numberOfFaces)
        error('mean_curvature_off:InvalidCounts', ...
              'Invalid vertex or face count in the OFF header.');
    end

    cursor = 5; % Token 4 is the number of edges (not needed here).
    requiredVertexTokens = 3 * numberOfVertices;
    if cursor + requiredVertexTokens - 1 > numel(tokens)
        error('mean_curvature_off:IncompleteVertices', ...
              'The OFF file ends before all vertices are defined.');
    end

    coordinates = str2double(tokens(cursor:cursor + requiredVertexTokens - 1));
    if any(~isfinite(coordinates))
        error('mean_curvature_off:InvalidVertices', ...
              'A vertex coordinate is not a finite number.');
    end
    V = reshape(coordinates, 3, numberOfVertices)';
    cursor = cursor + requiredVertexTokens;

    F = zeros(numberOfFaces, 3);
    for faceIndex = 1:numberOfFaces
        if cursor > numel(tokens)
            error('mean_curvature_off:IncompleteFaces', ...
                  'The OFF file ends before all faces are defined.');
        end
        verticesInFace = str2double(tokens{cursor});
        cursor = cursor + 1;
        if verticesInFace ~= 3
            error('mean_curvature_off:NonTriangularFace', ...
                  'Face %d has %g vertices; triangular faces are required.', ...
                  faceIndex, verticesInFace);
        end
        if cursor + 2 > numel(tokens)
            error('mean_curvature_off:IncompleteFace', ...
                  'Face %d is incomplete.', faceIndex);
        end
        F(faceIndex,:) = str2double(tokens(cursor:cursor + 2)) + 1;
        cursor = cursor + 3;
    end

    if any(~isfinite(F(:))) || any(F(:) ~= floor(F(:))) || ...
       any(F(:) < 1) || any(F(:) > numberOfVertices)
        error('mean_curvature_off:InvalidFaceIndices', ...
              'The OFF file contains an invalid face index.');
    end
end
