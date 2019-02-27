function phiGrad = cfdComputeGradientGauss0(phi,theMesh)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function calculates gradient based on green gauss method
%--------------------------------------------------------------------------

if(nargin<2)
    theMesh = cfdGetMesh;
end

theSize = size(phi);

theNumberOfComponents = theSize(2);
if(theNumberOfComponents > 3)
    echo('********* ERROR **********');
    exit;
end



%-----------------------------------------------------
% INTERIOR FACES contribution to gradient
%-----------------------------------------------------
iFaces = 1:theMesh.numberOfInteriorFaces;
iBFaces = theMesh.numberOfInteriorFaces+1:theMesh.numberOfFaces;

iElements = 1:theMesh.numberOfElements;
iBElements = theMesh.numberOfElements+1:theMesh.numberOfElements+theMesh.numberOfBFaces;

iOwners = [theMesh.faces(iFaces).iOwner]';
iNeighbours = [theMesh.faces(iFaces).iNeighbour]';

Sf = [theMesh.faces(iFaces).Sf]';
gf = [theMesh.faces(iFaces).gf]';
%-----------------------------------------------------
% Initialize phiGrad Array
%-----------------------------------------------------

phiGrad = zeros(theMesh.numberOfElements+theMesh.numberOfBElements,3,theNumberOfComponents);

for iComponent=1:theNumberOfComponents    
    phi_f = gf.*phi(iNeighbours,iComponent) + (1-gf).*phi(iOwners,iComponent);
    for iFace=iFaces
        phiGrad(iOwners(iFace),:,iComponent)     = phiGrad(iOwners(iFace),:,iComponent)     + phi_f(iFace)*Sf(iFace,:);
        phiGrad(iNeighbours(iFace),:,iComponent) = phiGrad(iNeighbours(iFace),:,iComponent) - phi_f(iFace)*Sf(iFace,:);
    end    
end

%-----------------------------------------------------
% BOUNDARY FACES contribution to gradient
%-----------------------------------------------------
iBOwners = [theMesh.faces(iBFaces).iOwner]';
phi_b = phi(iBElements,:);
Sb = [theMesh.faces(iBFaces).Sf]';
for iComponent=1:theNumberOfComponents
    %
    for k=1:theMesh.numberOfBFaces
        phiGrad(iBOwners(k),:,iComponent) = phiGrad(iBOwners(k),:,iComponent) + phi_b(k)*Sb(k,:);
    end
end


%-----------------------------------------------------
% Get Average Gradient by dividing with element volume
%-----------------------------------------------------
volumes = [theMesh.elements(iElements).volume]';
for iComponent=1:theNumberOfComponents
    for iElement =1:theMesh.numberOfElements
        phiGrad(iElement,:,iComponent) = phiGrad(iElement,:,iComponent)/volumes(iElement);
    end
end
%-----------------------------------------------------
% Set boundary Gradient equal to associated element
% Gradient
%-----------------------------------------------------
phiGrad(iBElements,:,:) = phiGrad(iBOwners,:,:);


end