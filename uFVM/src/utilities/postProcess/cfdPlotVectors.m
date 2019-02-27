function cfdPlotVectors(X, Y, Z, U, V, W, S)

q = quiver3(X, Y, Z, U, V, W, S);

%// Compute the magnitude of the vectors
mags = sqrt(sum(cat(2, q.UData(:), q.VData(:), ...
            reshape(q.WData, numel(q.UData), [])).^2, 2));

%// Get the current colormap
currentColormap = colormap(gca);

%// Now determine the color to make each arrow using a colormap
[~, ~, ind] = histcounts(mags, size(currentColormap, 1));

%// Now map this to a colormap to get RGB
cmap = uint8(ind2rgb(ind(:), currentColormap) * 255);
cmap(:,:,4) = 255;
cmap = permute(repmat(cmap, [1 3 1]), [2 1 3]);

%// We repeat each color 3 times (using 1:3 below) because each arrow has 3 vertices
set(q.Head, ...
    'ColorBinding', 'interpolated', ...
    'ColorData', reshape(cmap(1:3,:,:), [], 4).');   %'

%// We repeat each color 2 times (using 1:2 below) because each tail has 2 vertices
set(q.Tail, ...
    'ColorBinding', 'interpolated', ...
    'ColorData', reshape(cmap(1:2,:,:), [], 4).');