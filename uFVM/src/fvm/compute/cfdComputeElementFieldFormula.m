function elementField = cfdComputeElementFieldFormula(theFormula)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================

global time;
%
theMesh = cfdGetMesh;
%
%---------------------------------------------------
% update interior elements
%---------------------------------------------------
%
theNumberOfElements = theMesh.numberOfElements;
theNumberOfInteriorFaces = theMesh.numberOfInteriorFaces;
theNumberOfFaces = theMesh.numberOfFaces;

theNumberOfBElements =  theMesh.numberOfBElements;
theElementCentroids = [theMesh.elements(:).centroid]';
theBElementCentroids = [theMesh.faces(theNumberOfInteriorFaces+1:theNumberOfFaces).centroid]';

theCentroids = [theElementCentroids;theBElementCentroids];

x=theCentroids(:,1);
y=theCentroids(:,2);
z=theCentroids(:,3);

%
% Evaluate the Field
%
phi = eval(theFormula);
elementField = phi .* ones(theNumberOfElements+theNumberOfBElements,1);
%

