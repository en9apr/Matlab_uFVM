function theMeshField = cfdSetupMeshField(theUserName,theLocale,theType,theTimeStep)

%
%=================================================cfdSetupMeshField==
%  written by the CFD Group @ AUB, Fall 2006
%===================================================
%

if(nargin<4)
    theTimeStep = 'Step0'; % Step:0, Step:-1, Step:-2
end
if(nargin<3)
    theType = 'Scalar';
end
if(nargin<2)
    theLocale = 'Elements';
end
%
theMesh = cfdGetMesh;
%
% Name
%
theName = cfdConvertName(theUserName);
theMeshField.userName = theUserName;
theMeshField.name = theName;
%
% Properties
%
theMeshField.type = theType;
theMeshField.locale = theLocale;
%
% Locale
%
theInteriorArraySize=0;
theBoundaryArraySize=0;
if(strcmp(theLocale,'Elements'))
    theInteriorArraySize = theMesh.numberOfElements;
    theBoundaryArraySize = theMesh.numberOfBElements;
elseif(strcmp(theLocale,'Faces'))
    theInteriorArraySize = theMesh.numberOfInteriorFaces;
    theBoundaryArraySize = theMesh.numberOfBFaces;
elseif(strcmp(theLocale,'Nodes'))
    theInteriorArraySize = theMesh.numberOfNodes;
    theBoundaryArraySize = 0;
end
%
% Type
%
if(strcmp(theType,'Scalar'))
    theMeshField.phi = zeros(theInteriorArraySize+theBoundaryArraySize,1);
elseif(strcmp(theType,'Vector'))
     theMeshField.phi = zeros(theInteriorArraySize+theBoundaryArraySize,3);   
elseif(strcmp(theType,'Tensor'))
     theMeshField.phi = zeros(theInteriorArraySize+theBoundaryArraySize,9);      
end

%
% Defaults
%
cfdSetMeshField(theMeshField,theTimeStep);
%
