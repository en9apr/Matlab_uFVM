function cfdPlotPatchVelocity(theFluidName,iPatches,scale,iFigure)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================
%
theMesh = cfdGetMesh;

if(nargin==0)
    theFluid = cfdGetFluidUsingIndex(1);
    iPatches = 1:theMesh.numberOfPatches;
    iFigure = 1;
    scale = 0.5;
elseif(nargin==1)
    theFluid = cfdGetFluidUsingName(theFluidName);
    iPatches = 1:theMesh.numberOfPatches;
    iFigure = 1;
    scale = 0.5;
elseif(nargin==2)
    theFluid = cfdGetFluidUsingName(theFluidName);
    iFigure = 1;
    scale = 0.5;
elseif(nargin==3)
    theFluid = cfdGetFluidUsingName(theFluidName);
    iFigure = 1;
else
    theFluid = cfdGetFluidUsingName(theFluidName);
end

figure(iFigure);

%
theFluidTag = theFluid.tag;
theVelocityField = cfdGetMeshField(['Velocity' theFluidTag]);
theFacesVelocity = cfdInterpolateFromElementsToFaces('Average',theVelocityField.phi);
%
cfdPlotSurfaceMesh(iPatches,0,iFigure);
%

for iPatch=iPatches
    theBoundary = theMesh.boundaries(iPatch);
    numberOfBFaces = theBoundary.numberOfBFaces;
    %
    iFaceStart = theBoundary.startFace;
    iFaceEnd = iFaceStart+numberOfBFaces-1;
    iBFaces = iFaceStart:iFaceEnd;
    
    x = [];
    y = [];
    z = [];
    
    vx = [];
    vy = [];
    vz = [];
    
    for iBFace=iBFaces
        x = [x theMesh.faces(iBFace).centroid(1)];
        y = [y theMesh.faces(iBFace).centroid(2)];
        z = [z theMesh.faces(iBFace).centroid(3)];
        %
        vx = [vx theFacesVelocity(iBFace,1)];
        vy = [vy theFacesVelocity(iBFace,2)];
        vz = [vz theFacesVelocity(iBFace,3)];
    end
    
    hold on;
    quiver3(x,y,z,vx,vy,vz,scale);
    
end
%
axis equal;
colorbar;


