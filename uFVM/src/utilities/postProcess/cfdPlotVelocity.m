function cfdPlotVelocity(varargin)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================
%

% Default properties
faceAlpha = 0.01;
edgeAlpha = 0.1;
vectorSkip = 1;
vectorScale = 1;

for iProperty=1:2:length(varargin)
    if strcmp(varargin{iProperty},'faceAlpha')
        faceAlpha = varargin{iProperty+1};
    elseif strcmp(varargin{iProperty},'edgeAlpha')
        edgeAlpha = varargin{iProperty+1};
    elseif strcmp(varargin{iProperty},'vectorSkip')
        vectorSkip = varargin{iProperty+1};
    elseif strcmp(varargin{iProperty},'vectorScale')
        vectorScale = varargin{iProperty+1};
    end
end

figure;

axis off
axis equal

theMesh = cfdGetMesh;

theVelocityField = cfdGetMeshField('U');

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

cfdDrawMesh(v, f, 'FaceAlpha', faceAlpha, 'EdgeAlpha', edgeAlpha);

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

if vectorSkip==0
    for iElement=1:theNumberOfElements+theNumberOfBElements
        vx(iElement,1) = theVelocityField.phi(iElement,1);
        vy(iElement,1) = theVelocityField.phi(iElement,2);
        vz(iElement,1) = theVelocityField.phi(iElement,3);
    end
else
    vectorSkip = vectorSkip + 1;
    for iElement=1:vectorSkip:theNumberOfElements+theNumberOfBElements
        vx(iElement,1) = theVelocityField.phi(iElement,1);
        vy(iElement,1) = theVelocityField.phi(iElement,2);
        vz(iElement,1) = theVelocityField.phi(iElement,3);
    end
end
%

%
hold on;
colormap jet;
cfdPlotVectors(x,y,z,vx,vy,vz,vectorScale);
axis equal;

if min(cfdMagnitude([vx vy vz]))~=max(cfdMagnitude([vx vy vz]))
    caxis([min(cfdMagnitude([vx vy vz])) max(cfdMagnitude([vx vy vz]))]);
end

colorbar;

rotate3d;

set(gca,'Clipping','off');


