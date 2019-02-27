function cfdSetupVFFields

theNumberOfFluids = cfdGetNumberOfFluids;
for iFluid=1:theNumberOfFluids
    theVFName = ['VF' cfdGetFluidTagUsingIndex(iFluid)];
    theVFField = cfdGetModel(theVFName);
    if(isempty(theVFField))
        cfdSetupProperty(theVFName,'constant','1');
        %
        cfdSetupMeshField(theVFName,'Elements');
    end
    cfdInitializeProperty(theVFName);
    %
    theVFField = cfdGetMeshField(theVFName);
    VF_f = cfdSetupMeshField(theVFName,'Faces');
    phi_f = cfdInterpolateFromElementsToFaces('Average',theVFField.phi);
    VF_f.phi = phi_f;
    cfdSetMeshField(VF_f);

end
