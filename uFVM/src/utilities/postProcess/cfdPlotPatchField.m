function cfdPlotPatchField(theFieldName,varargin)

theMesh = cfdGetMesh;

if(numel(varargin)==0)
    iPatches = 1:theMesh.numberOfPatches;
    color = {'FaceAlpha', 1, 'EdgeAlpha', .1};    
elseif(numel(varargin)==1)
    iPatches = varargin{1};
    color = {'FaceAlpha', 1, 'EdgeAlpha', .1};
else
    iPatches = varargin{1};
    color = varargin{2};
end

theFigure = figure;
theFigure.Name = theFieldName;
colormap jet;

title(theFieldName);

theFluidName = cfdGetFluidName(theFieldName);

if strcmp(theFieldName,['MachNumber:' theFluidName])
    theMeshField = cfdGetMeshField(['Velocity:' theFluidName]);
    Vmag = cfdMagnitude(theMeshField.phi);
    
    theTemperatureField = cfdGetMeshField(['Temperature:' theFluidName]);
    T = theTemperatureField.phi;
    
    theGasConstantField = cfdGetMeshField(['GasConstant:' theFluidName]);
    R = theGasConstantField.phi;
    
    theSpecificHeatRatioField = cfdGetMeshField(['SpecificHeatRatio:' theFluidName]);
    k = theSpecificHeatRatioField.phi;
    
    theMeshField.phi = Vmag./sqrt(k.*R.*T);
    theMeshField.type = 'Scalar';
else
    theMeshField = cfdGetMeshField(theFieldName);
end

if(strcmp(theMeshField.type,'Vector'))    
    phiNodes = cfdInterpolateFromElementsToNodes(sqrt(dot(theMeshField.phi(:,:)',theMeshField.phi(:,:)')));
    cfdPlotSurfaceMesh(iPatches,phiNodes,theFigure);
    colorbar;    
else       
    phiNodes = cfdInterpolateFromElementsToNodes(theMeshField.phi(:,1));
    cfdPlotSurfaceMesh(iPatches,phiNodes,theFigure);
    colorbar;
end