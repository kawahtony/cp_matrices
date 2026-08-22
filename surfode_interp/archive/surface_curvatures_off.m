function [H, K, k1, k2, V, F, raw] = surface_curvatures_off( ...
    filename, smoothingIterations, smoothingStrength)
%SURFACE_CURVATURES_OFF Compute surface curvatures on a triangular OFF mesh.
%
%   [H,K,K1,K2,V,F] = SURFACE_CURVATURES_OFF('mesh.off') returns signed
%   mean curvature H, Gaussian curvature K, principal curvatures K1 >= K2,
%   vertices V, and triangular faces F.
%
%   [...] = SURFACE_CURVATURES_OFF('mesh.off', N, LAMBDA) applies N
%   positive-cotangent neighbor-smoothing passes with blend parameter
%   LAMBDA. Recommended values are N = 5 to 15 and LAMBDA = 0.5.
%   Set N = 0 to disable smoothing.
%
%   [H,K,K1,K2,V,F,RAW] also returns the unfiltered quantities in:
%       RAW.mean
%       RAW.gaussian
%       RAW.k1
%       RAW.k2
%       RAW.meanNormal
%
%   Conventions and formulas:
%   * H is positive on convex regions when face normals point outward.
%   * K is computed from the angle defect divided by barycentric area.
%   * K1,2 = H +/- sqrt(max(H^2-K,0)).
%   * K is independent of face orientation; H, K1, and K2 are not.
%   * Curvature units are 1/length for H,K1,K2 and 1/length^2 for K.
%
%   At a true mesh boundary, pi replaces 2*pi in the angle defect. That
%   boundary defect includes the discrete geodesic-curvature contribution.

    if nargin < 1
        error('surface_curvatures_off:MissingInput', ...
              'Supply the path to an OFF file.');
    end
    if nargin < 2
        smoothingIterations = 0;
    end
    if nargin < 3
        smoothingStrength = 0.5;
    end
    validateattributes(smoothingIterations, {'numeric'}, ...
                       {'scalar', 'real', 'finite', 'integer', 'nonnegative'}, ...
                       mfilename, 'smoothingIterations', 2);
    validateattributes(smoothingStrength, {'numeric'}, ...
                       {'scalar', 'real', 'finite', '>=', 0, '<=', 1}, ...
                       mfilename, 'smoothingStrength', 3);

    [V, F] = read_triangular_off(filename);
    numberOfVertices = size(V, 1);
    numberOfFaces = size(F, 1);

    p1 = V(F(:,1),:);
    p2 = V(F(:,2),:);
    p3 = V(F(:,3),:);
    faceNormal = cross(p2 - p1, p3 - p1, 2);
    doubleArea = sqrt(sum(faceNormal.^2, 2));

    meshScale = norm(max(V, [], 1) - min(V, [], 1));
    areaTolerance = max(realmin, 1e-14 * meshScale^2);
    nondegenerate = doubleArea > areaTolerance;
    if any(~nondegenerate)
        warning('surface_curvatures_off:DegenerateFaces', ...
                '%d degenerate face(s) were ignored.', ...
                sum(~nondegenerate));
    end

    % Corner cotangents and angles. atan2 is stable near 0 and pi.
    cot1 = zeros(numberOfFaces, 1);
    cot2 = zeros(numberOfFaces, 1);
    cot3 = zeros(numberOfFaces, 1);
    angle1 = zeros(numberOfFaces, 1);
    angle2 = zeros(numberOfFaces, 1);
    angle3 = zeros(numberOfFaces, 1);

    dot1 = row_dot(p2 - p1, p3 - p1);
    dot2 = row_dot(p1 - p2, p3 - p2);
    dot3 = row_dot(p1 - p3, p2 - p3);
    cot1(nondegenerate) = dot1(nondegenerate) ./ doubleArea(nondegenerate);
    cot2(nondegenerate) = dot2(nondegenerate) ./ doubleArea(nondegenerate);
    cot3(nondegenerate) = dot3(nondegenerate) ./ doubleArea(nondegenerate);
    angle1(nondegenerate) = atan2(doubleArea(nondegenerate), dot1(nondegenerate));
    angle2(nondegenerate) = atan2(doubleArea(nondegenerate), dot2(nondegenerate));
    angle3(nondegenerate) = atan2(doubleArea(nondegenerate), dot3(nondegenerate));

    % Cotangent matrix and stiffness matrix for the mean-curvature normal.
    i = [F(:,2); F(:,3); F(:,3); F(:,1); F(:,1); F(:,2)];
    j = [F(:,3); F(:,2); F(:,1); F(:,3); F(:,2); F(:,1)];
    weights = [cot1; cot1; cot2; cot2; cot3; cot3];
    W = sparse(i, j, weights, numberOfVertices, numberOfVertices);
    L = spdiags(sum(W, 2), 0, numberOfVertices, numberOfVertices) - W;

    faceArea = 0.5 * doubleArea;
    vertexArea = accumarray(F(:), repmat(faceArea / 3, 3, 1), ...
                            [numberOfVertices, 1], @sum, 0);
    validVertex = vertexArea > areaTolerance;

    meanNormal = nan(numberOfVertices, 3);
    meanNormal(validVertex,:) = bsxfun( ...
        @rdivide, L(validVertex,:) * V, 4 * vertexArea(validVertex));

    % Area-weighted vertex normals.
    vertexNormal = zeros(numberOfVertices, 3);
    repeatedFaceNormal = repmat(faceNormal, 3, 1);
    for coordinate = 1:3
        vertexNormal(:,coordinate) = accumarray( ...
            F(:), repeatedFaceNormal(:,coordinate), ...
            [numberOfVertices, 1], @sum, 0);
    end
    normalLength = sqrt(sum(vertexNormal.^2, 2));
    hasNormal = normalLength > 0;
    vertexNormal(hasNormal,:) = bsxfun( ...
        @rdivide, vertexNormal(hasNormal,:), normalLength(hasNormal));

    % Detect true boundary and nonmanifold edges without forming a full
    % vertex-by-vertex matrix.
    edgeStart = [F(:,1); F(:,2); F(:,3)];
    edgeEnd = [F(:,2); F(:,3); F(:,1)];
    lowerVertex = min(edgeStart, edgeEnd);
    upperVertex = max(edgeStart, edgeEnd);
    edgeIncidence = sparse(lowerVertex, upperVertex, ...
                           ones(size(lowerVertex)), ...
                           numberOfVertices, numberOfVertices);
    [edgeI, edgeJ, incidenceCount] = find(edgeIncidence);
    isBoundary = false(numberOfVertices, 1);
    boundaryEdge = incidenceCount == 1;
    isBoundary([edgeI(boundaryEdge); edgeJ(boundaryEdge)]) = true;
    if any(incidenceCount > 2)
        warning('surface_curvatures_off:NonmanifoldEdges', ...
                '%d edge(s) have more than two incident faces.', ...
                sum(incidenceCount > 2));
    end
    clear edgeIncidence edgeI edgeJ incidenceCount

    % Angle-defect Gaussian curvature.
    angleSum = accumarray(F(:), [angle1; angle2; angle3], ...
                          [numberOfVertices, 1], @sum, 0);
    targetAngle = 2 * pi * ones(numberOfVertices, 1);
    targetAngle(isBoundary) = pi;
    gaussianRaw = nan(numberOfVertices, 1);
    gaussianRaw(validVertex) = ...
        (targetAngle(validVertex) - angleSum(validVertex)) ...
        ./ vertexArea(validVertex);

    meanRaw = sum(meanNormal .* vertexNormal, 2);
    meanRaw(~hasNormal) = NaN;
    [k1Raw, k2Raw] = principal_from_mean_gaussian(meanRaw, gaussianRaw);

    raw = struct('mean', meanRaw, ...
                 'gaussian', gaussianRaw, ...
                 'k1', k1Raw, ...
                 'k2', k2Raw, ...
                 'meanNormal', meanNormal);

    K = gaussianRaw;
    if smoothingIterations > 0 && smoothingStrength > 0
        positiveW = spfun(@(weight) max(weight, 0), W);
        meanNormal = smooth_cotangent_field( ...
            meanNormal, positiveW, smoothingIterations, smoothingStrength);
        K = smooth_cotangent_field( ...
            K, positiveW, smoothingIterations, smoothingStrength);

        % Preserve the area-weighted Gaussian-curvature integral.
        finiteRaw = validVertex & isfinite(gaussianRaw);
        finiteFiltered = validVertex & isfinite(K);
        integralBefore = sum(vertexArea(finiteRaw) .* gaussianRaw(finiteRaw));
        integralAfter = sum(vertexArea(finiteFiltered) .* K(finiteFiltered));
        filteredArea = sum(vertexArea(finiteFiltered));
        if filteredArea > 0
            K(finiteFiltered) = K(finiteFiltered) + ...
                (integralBefore - integralAfter) / filteredArea;
        end
    end

    H = sum(meanNormal .* vertexNormal, 2);
    H(~hasNormal) = NaN;
    K(~validVertex) = NaN;
    [k1, k2] = principal_from_mean_gaussian(H, K);
end


function field = smooth_cotangent_field(field, positiveW, iterations, lambda)
%SMOOTH_COTANGENT_FIELD Stable filtering of a scalar or vector vertex field.
    finiteVertex = all(isfinite(field), 2);

    for iteration = 1:iterations
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


function [k1, k2] = principal_from_mean_gaussian(H, K)
%PRINCIPAL_FROM_MEAN_GAUSSIAN Recover ordered principal curvatures.
%   A negative discriminant can occur because H and K are independent
%   discrete estimates. Clamping it to zero gives the closest real pair.
    discriminant = H.^2 - K;
    invalid = ~isfinite(discriminant);
    discriminant = max(discriminant, 0);
    discriminant(invalid) = NaN;
    halfDifference = sqrt(discriminant);
    k1 = H + halfDifference;
    k2 = H - halfDifference;
end


function d = row_dot(a, b)
%ROW_DOT Dot product of corresponding matrix rows.
    d = sum(a .* b, 2);
end


function [V, F] = read_triangular_off(filename)
%READ_TRIANGULAR_OFF Read vertices and triangular faces from an ASCII OFF file.
    fileId = fopen(filename, 'r');
    if fileId == -1
        error('surface_curvatures_off:CannotOpenFile', ...
              'Could not open OFF file: %s', filename);
    end
    closeFile = onCleanup(@() fclose(fileId));

    contents = fread(fileId, '*char')';
    contents = regexprep(contents, '#[^\r\n]*', ' ');
    tokens = regexp(contents, '\S+', 'match');

    if isempty(tokens) || ~strcmpi(tokens{1}, 'OFF')
        error('surface_curvatures_off:InvalidHeader', ...
              'The file must begin with the ASCII OFF header.');
    end
    if numel(tokens) < 4
        error('surface_curvatures_off:IncompleteFile', ...
              'The OFF file does not contain valid mesh counts.');
    end

    numberOfVertices = str2double(tokens{2});
    numberOfFaces = str2double(tokens{3});
    if ~isfinite(numberOfVertices) || ~isfinite(numberOfFaces) || ...
       numberOfVertices < 1 || numberOfFaces < 1 || ...
       numberOfVertices ~= floor(numberOfVertices) || ...
       numberOfFaces ~= floor(numberOfFaces)
        error('surface_curvatures_off:InvalidCounts', ...
              'Invalid vertex or face count in the OFF header.');
    end

    cursor = 5;
    requiredVertexTokens = 3 * numberOfVertices;
    if cursor + requiredVertexTokens - 1 > numel(tokens)
        error('surface_curvatures_off:IncompleteVertices', ...
              'The OFF file ends before all vertices are defined.');
    end

    coordinates = str2double(tokens(cursor:cursor + requiredVertexTokens - 1));
    if any(~isfinite(coordinates))
        error('surface_curvatures_off:InvalidVertices', ...
              'A vertex coordinate is not a finite number.');
    end
    V = reshape(coordinates, 3, numberOfVertices)';
    cursor = cursor + requiredVertexTokens;

    F = zeros(numberOfFaces, 3);
    for faceIndex = 1:numberOfFaces
        if cursor > numel(tokens)
            error('surface_curvatures_off:IncompleteFaces', ...
                  'The OFF file ends before all faces are defined.');
        end
        verticesInFace = str2double(tokens{cursor});
        cursor = cursor + 1;
        if verticesInFace ~= 3
            error('surface_curvatures_off:NonTriangularFace', ...
                  'Face %d has %g vertices; triangular faces are required.', ...
                  faceIndex, verticesInFace);
        end
        if cursor + 2 > numel(tokens)
            error('surface_curvatures_off:IncompleteFace', ...
                  'Face %d is incomplete.', faceIndex);
        end
        F(faceIndex,:) = str2double(tokens(cursor:cursor + 2)) + 1;
        cursor = cursor + 3;
    end

    if any(~isfinite(F(:))) || any(F(:) ~= floor(F(:))) || ...
       any(F(:) < 1) || any(F(:) > numberOfVertices)
        error('surface_curvatures_off:InvalidFaceIndices', ...
              'The OFF file contains an invalid face index.');
    end
end
