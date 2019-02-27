function centroids = cfdGetElementCentroids;
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================

m=cfdGetMesh;
centroids = fvmCellsToArray({m.elements.centroid});