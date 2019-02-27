function cfdPlotStreamlines(theFluidName,skip,iFigure)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================
%
if(nargin == 0)
    theFluid = cfdGetFluidUsingIndex(1);
    iFigure = 1;
    skip = 1;
elseif(nargin == 1)
    theFluid = cfdGetFluidUsingName(theFluidName);
    skip = 1;
    iFigure = 1;
elseif(nargin == 2)
    theFluid= cfdGetFluidUsingName(theFluidName);
    iFigure = 1;
elseif(nargin == 3)
    theFluid= cfdGetFluidUsingName(theFluidName);
end

figure(iFigure)

axis off
axis equal

theMesh = cfdGetMesh;
%
theFluidTag = theFluid.tag;

theVelocityField = cfdGetMeshField(['Velocity' theFluidTag]);
%

theNumberOfElements = theMesh.numberOfElements;
theNumberOfBElements = theMesh.numberOfBElements;
theNumberOfNodes = theMesh.numberOfNodes;
theNumberOfFaces = theMesh.numberOfFaces;

for iNode=1:theNumberOfNodes
    v(iNode,:) = theMesh.nodes(iNode).centroid';
end

for iFace=1:theNumberOfFaces
    f{iFace,:} = theMesh.faces(iFace).iNodes;
end

f = cfdConvertCelltoMatrix(f);

cfdDrawMesh(v, f, 'FaceAlpha', 0.1, 'EdgeAlpha', .4);

view(30,40);
%
theNumberOfInteriorFaces = theMesh.numberOfInteriorFaces;

theElementCentroids = [theMesh.elements(:).centroid]';
theBElementCentroids = [theMesh.faces(theNumberOfInteriorFaces+1:theNumberOfFaces).centroid]';

theCentroids = [theElementCentroids;theBElementCentroids];

x = theCentroids(:,1);
y = theCentroids(:,2);
z = theCentroids(:,3);

%

vx = zeros(theNumberOfElements+theNumberOfBElements,1);
vy = zeros(theNumberOfElements+theNumberOfBElements,1);
vz = zeros(theNumberOfElements+theNumberOfBElements,1);

for iElement=1:skip:theNumberOfElements+theNumberOfBElements    
    vx(iElement,1) = theVelocityField.phi(iElement,1);
    vy(iElement,1) = theVelocityField.phi(iElement,2);
    vz(iElement,1) = theVelocityField.phi(iElement,3);    
end
%

%
hold on;
XYZ = stream3(x,y,z,vx,vy,vz,0,0,-0.03);
streamline(x,y,z,vx,vy,vz,0,0,-0.03);
axis equal;
colorbar;

rotate3d;


