function effDiv = cfdComputeEffectiveDivergence
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017 
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function computes the effective divergence of the continuity
%   equation
%--------------------------------------------------------------------------

% get mesh attributes
theMesh = cfdGetMesh;
theNumberOfElements = theMesh.numberOfElements;

% Get mdot_f field
theMdotField = cfdGetMeshField('mdot_f', 'Faces');
mdot_f = theMdotField.phi;

% Initialize effective divergence
effDiv = zeros(theNumberOfElements,1);

% Interior Faces Contribution
theNumberOfInteriorFaces = cfdGetNumberOfInteriorFaces;
iFaces = 1:theNumberOfInteriorFaces;
owners = [theMesh.faces(iFaces).iOwner]';
neighbours = [theMesh.faces(iFaces).iNeighbour]';
for iFace=1:theNumberOfInteriorFaces
    iOwner = owners(iFace);
    iNeighbour = neighbours(iFace);
    %
    effDiv(iOwner)     = effDiv(iOwner)     + mdot_f(iFace);
    effDiv(iNeighbour) = effDiv(iNeighbour) - mdot_f(iFace);
end

% Boundary Faces Contribution
theNumberOfPatches = cfdGetNumberOfPatches;
for iPatch=1:theNumberOfPatches    
    theBoundary = theMesh.boundaries(iPatch);
    numberOfBFaces = theBoundary.numberOfBFaces;
    
    % cfdGetBoundaryIndex
    iFaceStart = theBoundary.startFace;
    iFaceEnd = iFaceStart+numberOfBFaces-1;
    iBFaces = iFaceStart:iFaceEnd;

    owners = [theMesh.faces(iBFaces).iOwner]';
    mdot_b = theMdotField.phi(iBFaces);
    
    for iBFace=1:numberOfBFaces
        iOwner = owners(iBFace);
        effDiv(iOwner) = effDiv(iOwner) + mdot_b(iBFace);
    end    
end

end
