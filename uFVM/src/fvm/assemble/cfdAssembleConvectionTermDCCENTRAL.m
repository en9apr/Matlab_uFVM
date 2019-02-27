function theFluxes = cfdAssembleConvectionTermDCCENTRAL(theEquationName,theTerm,theFluxes,iComponent)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================
if(nargin == 3)
    iComponent = 1;
end

theMesh = cfdGetMesh;

numberOfInteriorFaces = theMesh.numberOfInteriorFaces;
iFaces = 1:numberOfInteriorFaces;

theEquationMeshField = cfdGetMeshField(theEquationName);
phi = theEquationMeshField.phi;
phiGrad = theEquationMeshField.phiGradient(:,:,1);

% Get the term coefficient
if strcmp(theEquationName, 'U')
    thePsiField = cfdGetMeshField('mdot_f', 'Faces');
else
    thePsiField = cfdGetMeshField(['psi_',theEquationName,'eq'], 'Faces');
end
psi = thePsiField.phi(iFaces);

iOwners = [theMesh.faces(iFaces).iOwner]';
iNeighbours = [theMesh.faces(iFaces).iNeighbour]';
pos = zeros(size(psi));
pos(psi>0) = 1;

iUpwind = pos.*iOwners + (1-pos).*iNeighbours;

phiGradf = cfdInterpolateGradientsFromElementsToInteriorFaces('Average:Corrected',phiGrad,phi);

rC = [theMesh.elements(iUpwind).centroid]';
rf = [theMesh.faces(iFaces).centroid]';

rCf = rf - rC;

corr = psi .* dot(phiGradf',rCf')';

theFluxes.FLUXTf(iFaces) = theFluxes.FLUXTf(iFaces) +  corr;
