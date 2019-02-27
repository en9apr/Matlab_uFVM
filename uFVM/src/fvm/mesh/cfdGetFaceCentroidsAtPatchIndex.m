function centroids = cfdGetFacesCentroidsAtPatchIndex(thePatchIndex)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================

m=cfdGetMesh;
centroids = fvmCellsToArray({m.faces(find([m.faces(:).patchIndex]==thePatchIndex)).centroid});