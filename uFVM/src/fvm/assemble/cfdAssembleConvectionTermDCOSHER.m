function theFluxes = cfdAssembleConvectionTermDCOSHER(theEquationName,theTerm,theFluxes)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================

theMesh = cfdGetMesh;
theEquationMeshField = cfdGetMeshField(theEquationName);
phi = theEquationMeshField.phi;
phiGrad = theEquationMeshField.phiGradient;
theFluidTag = cfdGetFluidTag(theEquationName);
theMdotName = ['Mdot' theFluidTag];
theMdotField = cfdGetMeshField(theMdotName,'Faces');


numberOfInteriorFaces = theMesh.numberOfInteriorFaces;
for iFaces = 1:numberOfInteriorFaces;



    mdot_f = theMdotField.phi(iFaces);
    iOwners = [theMesh.faces(iFaces).iOwner];
    iNeighbours = [theMesh.faces(iFaces).iNeighbour];
    pos = zeros(size(mdot_f));
    pos(find(mdot_f>0))=1;

    iUpwind = pos.*iOwners + (1-pos).*iNeighbours;
    iDownwind = (1-pos).*iOwners + pos.*iNeighbours;

    %get the upwind gradient at the interior faces
    phiGradUpwindf = phiGrad(:,iUpwind);
    %interpolated gradient to interior faces
    phiGradf_all = cfdInterpolateGradientsFromElementsToInteriorFaces(phiGrad,phi);
    phiGradf=phiGradf_all(:,iFaces);

    phiC=phi(iUpwind);
    phiD=phi(iDownwind);
    phiGradC=phiGrad(:,iUpwind);
    rC = [theMesh.elements(iUpwind).centroid];
    rD = [theMesh.elements(iDownwind).centroid];
    rCD = rD-rC;
    phiU= phiD -dot(2*phiGradC,rCD);
    if(abs(phiD-phiU)>1e-9)
        phiNormC = (phiC-phiU)/(phiD-phiU);
    else
        phiNormC = 0;
    end
    rf = [theMesh.faces(iFaces).centroid];
    rCf = rf-rC;

     if((phiNormC>=0) && (phiNormC<=2/3)) %SOU
        corr = mdot_f * (dot(2*phiGradUpwindf - phiGradf,rCf)-phiC);
        theFluxes.FLUXTf(iFaces) = theFluxes.FLUXTf(iFaces) - corr;
     elseif((phiNormC>2/3) && (phiNormC<=1)) %DW
        fluxC1f = mdot_f.*(1-pos);
        fluxC2f = mdot_f.*pos;
        tFLUXVf = 0;
        phiMdotSOU = fluxC1f* phi(iOwners) + fluxC2f * phi(iNeighbours) + FLUXV
        theFluxes.FLUXTf(iFaces) = -(phiMdotSOU-mdot_f*phiC) ;
     else %UW
        fluxC1f = mdot_f.*pos;
        fluxC2f = mdot_f.*(1-pos);
        theFluxes.FLUXVf(iFaces) = 0;
        phiMdotSOU = fluxC1f* phi(iOwners) + fluxC2f * phi(iNeighbours) + FLUXV
        theFluxes.FLUXTf(iFaces) = -(phiMdotSOU-mdot_f*phiC) ;
     end

end
