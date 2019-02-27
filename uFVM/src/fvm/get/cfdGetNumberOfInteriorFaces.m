function theNumberOfInteriorFaces = cfdGetNumberOfInteriorFaces
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================

 
 m = cfdGetMesh;
 theNumberOfInteriorFaces = m.numberOfInteriorFaces;