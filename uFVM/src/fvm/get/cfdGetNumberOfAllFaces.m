function theNumberOfAllFaces = cfdGetNumberOfAllFaces
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================

m = cfdGetMesh;
theNumberOfAllFaces = length(m.faces);