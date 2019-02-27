function drhodp = cfdGetDrhodp(theFluidTag,varargin)

theMesh = cfdGetMesh;
numberOfElements = theMesh.numberOfElements;
numberOfBElements = theMesh.numberOfBElements;
numberOfFaces = theMesh.numberOfFaces;

if nargin==1
    theLocale = 'Elements';
    iLocale = 1:numberOfElements+numberOfBElements;
elseif nargin>1
    theLocale = varargin{1};
    if strcmp(theLocale,'Elements')
        iLocale = 1:numberOfElements+numberOfBElements;
    elseif strcmp(theLocale,'Faces')
        iLocale = 1:numberOfFaces;
    end
end
if nargin>2
    iLocale = varargin{2};
end

theDensityModel = cfdGetModel(['Density' theFluidTag]);
model = theDensityModel.model;

if strcmp(model,'Ideal Gas')    
    temperature = cfdGetMeshField(['Temperature' theFluidTag]);
    T = temperature.phi;
    gasConstant = cfdGetMeshField(['GasConstant' theFluidTag]);
    R = gasConstant.phi;
        
    drhodp = 1./(R.*T);
    
    if strcmp(theLocale,'Faces')
        drhodp = cfdInterpolateFromElementsToFaces('Average',drhodp);
    end
    drhodp = drhodp(iLocale);
elseif strcmp(model,'Real Gas')
    temperature = cfdGetMeshField(['Temperature' theFluidTag]);
    T = temperature.phi;
    gasConstant = cfdGetMeshField(['GasConstant' theFluidTag]);
    R = gasConstant.phi;
    
    drhodp = 1./(R.*T);
    if strcmp(theLocale,'Faces')
        drhodp = cfdInterpolateFromElementsToFaces('Average',drhodp);
    end
    drhodp = drhodp(iLocale);
end
