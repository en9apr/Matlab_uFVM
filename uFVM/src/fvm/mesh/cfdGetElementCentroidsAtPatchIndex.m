function centroids = cfdGetElementCentroidsAtPatchIndex(thePatchIndex);
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================

m=cfdGetMesh;
centroids = fvmCellsToArray({(m.elements([m.faces(find([m.faces(:).patchIndex]==thePatchIndex)).element1]).centroid)});