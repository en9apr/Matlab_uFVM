function mesh = cfdReadPolyMesh
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017 
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function reads poly mesh files from "constant/polyMesh" directory
%   and stores them in the database ~domain
%--------------------------------------------------------------------------

fprintf('\nReading Poly Mesh ...\n');

projectDirectory = pwd;
projectDirectory = strrep(projectDirectory, '\', '/');

pointsFile = [projectDirectory, '/constant/polyMesh/points'];
facesFile = [projectDirectory, '/constant/polyMesh/faces'];
ownerFile = [projectDirectory, '/constant/polyMesh/owner'];
neighbourFile = [projectDirectory, '/constant/polyMesh/neighbour'];
boundaryFile = [projectDirectory, '/constant/polyMesh/boundary'];
%
% Read Points
%
fpid = fopen(pointsFile, 'r');
for i=1:18
    fgetl(fpid);
end
numberOfPoints = fscanf(fpid, '%d',1);
fgetl(fpid);
fgetl(fpid);
for n = 1:numberOfPoints
    fscanf(fpid, '%c',1);
    x=fscanf(fpid, '%f',1);
    y=fscanf(fpid, '%f',1);
    z=fscanf(fpid, '%f',1);
    fscanf(fpid, '%c',1);
    fscanf(fpid, '%c',1);
    nodes(n).centroid = [x y z]';
    nodes(n).index = n;
    nodes(n).iFaces=[];
    nodes(n).iElements=[];
end
mesh.nodes = nodes;
mesh.numberOfNodes = numberOfPoints;

mesh.caseDirectory = 'constant/polyMesh';


%
%read Faces
%
ffid = fopen(facesFile, 'r');
%header = textscan(fpid,'%s',18,'EndOfLine','\n') 
numberOfFaces = textscan(ffid,'%d',1,'HeaderLines',18);
numberOfFaces = numberOfFaces{1};
header = fscanf(ffid, '%s', 1);

% read each set of measurements
for n = 1:numberOfFaces
    theNumberOfPoints = fscanf(ffid, '%d', 1);
    header = fscanf(ffid, '%c', 1);
    faces(n).iNodes=fscanf(ffid, '%d', theNumberOfPoints)'+1;
    header = fscanf(ffid, '%c', 1);
    faces(n).index=n;
    faces(n).iOwner = -1;
    faces(n).iNeighbour = -1;
end
mesh.numberOfFaces = numberOfFaces;


%
%read owners
%
 foid = fopen(ownerFile, 'r');
%header = textscan(fpid,'%s',18,'EndOfLine','\n') 
numberOfOwners = textscan(foid,'%d',1,'HeaderLines',18);
numberOfOwners = numberOfOwners{1};
header = fscanf(foid, '%s', 1);
for n=1:numberOfOwners
    faces(n).iOwner = fscanf(foid,'%d',1)+1;
    faces(n).iNeighbour=-1;
end
mesh.numberOfElements = max([faces.iOwner]);

%
%read neighbours
%
fnid = fopen(neighbourFile, 'r');
%header = textscan(fpid,'%s',18,'EndOfLine','\n') 
numberOfNeighbours = textscan(fnid,'%d',1,'HeaderLines',18);
numberOfNeighbours = numberOfNeighbours{1};
header = fscanf(fnid, '%s', 1);
for n=1:numberOfNeighbours
    faces(n).iNeighbour = fscanf(fnid,'%d',1)+1;
end

mesh.faces = faces;
mesh.numberOfInteriorFaces = numberOfNeighbours;

numberOfElements = max([faces.iNeighbour]);
mesh.numberOfElements = max(mesh.numberOfElements,numberOfElements);

%
%read boundaries
%
fbid = fopen(boundaryFile, 'r');
numberOfBoundaries = textscan(fbid,'%d',1,'HeaderLines',17);
numberOfBoundaries = numberOfBoundaries{1};
header = fscanf(fbid, '%c', 2);
for n=1:numberOfBoundaries
       boundaries(n).userName = fscanf(fbid, '%s', 1);
        boundaries(n).index=n;
        token='{';
     while(strcmp(token,'}')==false)
        token = fscanf(fbid, '%s', 1);
        if(strcmp(token,'type'))
            theType = fscanf(fbid, '%s', 1);
            boundaries(n).type = theType(1:length(theType)-1);   
        end
        if(strcmp(token,'startFace'))
              boundaries(n).startFace = fscanf(fbid, '%d', 1)+1;   
              
        end
        if(strcmp(token,'nFaces'))
              boundaries(n).numberOfBFaces = fscanf(fbid, '%d', 1);   
        end
        fgetl(fbid);
    end
end
mesh.boundaries = boundaries;
mesh.numberOfBoundaries = numberOfBoundaries;
mesh.numberOfPatches = numberOfBoundaries;


%
% construct elements
%
elements= [];
numberOfElements = mesh.numberOfElements;
for iElement=1:numberOfElements
   elements(iElement).index = iElement;
   elements(iElement).iNeighbours = [];
   elements(iElement).iFaces=[];
   elements(iElement).iNodes=[];
   elements(iElement).volume = 0.;
   elements(iElement).faceSign=[];

end

numberOfInteriorFaces = mesh.numberOfInteriorFaces;
for iFace=1:numberOfInteriorFaces
   iOwner = faces(iFace).iOwner;
   iNeighbour = faces(iFace).iNeighbour;
   %
   elements(iOwner).iFaces = [elements(iOwner).iFaces iFace];
   elements(iOwner).faceSign=[elements(iOwner).faceSign 1];
   elements(iOwner).iNeighbours = [elements(iOwner).iNeighbours iNeighbour];
   %
   elements(iNeighbour).iFaces = [elements(iNeighbour).iFaces iFace];
   elements(iNeighbour).faceSign=[elements(iNeighbour).faceSign -1];
   elements(iNeighbour).iNeighbours = [elements(iNeighbour).iNeighbours iOwner];
end

numberOfFaces = mesh.numberOfFaces;
for iBFace=numberOfInteriorFaces+1:numberOfFaces
   iOwner = mesh.faces(iBFace).iOwner;
   %
   elements(iOwner).iFaces = [elements(iOwner).iFaces iBFace];
   elements(iOwner).faceSign=[elements(iOwner).faceSign 1];
end

numberOfElements = mesh.numberOfElements;
for iElement=1:numberOfElements
   elements(iElement).numberOfNeighbours = length(elements(iElement).iNeighbours);
end

mesh.elements=elements;

mesh.numberOfBElements = numberOfFaces - numberOfInteriorFaces;
mesh.numberOfBFaces = numberOfFaces - numberOfInteriorFaces;

%
% setup Node connectivities
%
for iFace=1:numberOfFaces    
    iNodes = mesh.faces(iFace).iNodes;
    for iNode=iNodes
       mesh.nodes(iNode).iFaces = [mesh.nodes(iNode).iFaces iFace];
    end    
end

%
% setup Node connectivities
%
for iElement=1:numberOfElements
    
    iFaces = mesh.elements(iElement).iFaces;
    for iFace=iFaces
       iNodes = mesh.faces(iFace).iNodes;
       for iNode=iNodes
           if(sum(mesh.elements(iElement).iNodes == iNode)==0)
              mesh.elements(iElement).iNodes = [mesh.elements(iElement).iNodes iNode];
              mesh.nodes(iNode).iElements = [mesh.nodes(iNode).iElements iElement];
           end
       end
    end
    
end

theMesh = cfdProcessOpenFoamMesh(mesh);

% Store
cfdSetMesh(theMesh);





