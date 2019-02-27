function cfdAssembleIntoGlobalMatrixFaceFluxes(theFluxes)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017 
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function assembles algebraic equation coefficients from the
%   contribution of the face fluxes of the current term of the equation
%--------------------------------------------------------------------------


theMesh = cfdGetMesh;
numberOfInteriorFaces = theMesh.numberOfInteriorFaces;
numberOfFaces = theMesh.numberOfFaces;

theCoefficients = cfdGetCoefficients;

ac = theCoefficients.ac;
anb = theCoefficients.anb;
bc = theCoefficients.bc;
%
% Assemble fluxes of interior faces
%
for iFace = 1:numberOfInteriorFaces
    theFace                 = theMesh.faces(iFace);
    iOwner                  = theFace.iOwner;
    iOwnerNeighbourCoef     = theFace.iOwnerNeighbourCoef;
    iNeighbour              = theFace.iNeighbour;
    iNeighbourOwnerCoef     = theFace.iNeighbourOwnerCoef;
    % 
    %  assemble fluxes for owner cell
    %
    ac(iOwner)                       = ac(iOwner)                       + theFluxes.FLUXC1f(iFace);
    anb{iOwner}(iOwnerNeighbourCoef) = anb{iOwner}(iOwnerNeighbourCoef) + theFluxes.FLUXC2f(iFace);
    bc(iOwner)                       = bc(iOwner)                       - theFluxes.FLUXTf(iFace);
    % 
    %  assemble fluxes for neighbour cell
    %
    ac(iNeighbour)                       = ac(iNeighbour)                       - theFluxes.FLUXC2f(iFace);
    anb{iNeighbour}(iNeighbourOwnerCoef) = anb{iNeighbour}(iNeighbourOwnerCoef) - theFluxes.FLUXC1f(iFace);
    bc(iNeighbour)                       = bc(iNeighbour)                       + theFluxes.FLUXTf(iFace);
end
%
% assemble fluxes of boundary faces
%
for iBFace=numberOfInteriorFaces+1:numberOfFaces
    theBFace = theMesh.faces(iBFace);
    iOwner         = theBFace.iOwner;
    %
    %  assemble fluxes for owner cell
    %
    ac(iOwner) = ac(iOwner) + theFluxes.FLUXC1f(iBFace);
    bc(iOwner) = bc(iOwner) - theFluxes.FLUXTf(iBFace);
end

theCoefficients.ac = ac;
theCoefficients.anb = anb;
theCoefficients.bc = bc;

cfdSetCoefficients(theCoefficients);

end
