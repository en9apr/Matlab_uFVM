function effDiv = cfdComputeEffectiveDivergenceWithFluidTag(theTerm)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017 
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function claculates effective divergence
%--------------------------------------------------------------------------

theMesh = cfdGetMesh;

% Get the term coefficient
if strcmp(theTerm.variableName, 'U')
    thePsiField = cfdGetMeshField('mdot_f', 'Faces');
else
    thePsiField = cfdGetMeshField(['psi_',theTerm.variableName,'eq'], 'Faces');
end
psi = thePsiField.phi;

effDiv = zeros(theMesh.numberOfElements+theMesh.numberOfBFaces,1);
%
% Interior Faces Contribution
%
iFaces = 1:theMesh.numberOfInteriorFaces;
iOwners(iFaces,1) = [theMesh.faces(iFaces).iOwner]';
iNeighbours(iFaces,1) = [theMesh.faces(iFaces).iNeighbour]';

for iFace=iFaces
    iOwner = iOwners(iFace);
    iNeighbour = iNeighbours(iFace);
    
    effDiv(iOwner) = effDiv(iOwner) + psi(iFace);
    effDiv(iNeighbour) = effDiv(iNeighbour) - psi(iFace);
end

%
% Boundary Faces Contribution
%
iBFaces = theMesh.numberOfInteriorFaces+1:theMesh.numberOfFaces;
iOwners(iBFaces) = [theMesh.faces(iBFaces).iOwner]';

for iFace=iBFaces
    iOwner = iOwners(iFace);    
    effDiv(iOwner) = effDiv(iOwner) + psi(iFace);
end

end
