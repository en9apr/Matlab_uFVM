function sumVF = cfdSumVFFields
%
%
%
VF=cfdGetMeshField('VF_fluid01');

sumVF=zeros(size(VF.phi));

theNumberOfFluids = cfdGetNumberOfFluids;
for iFluid=1:theNumberOfFluids
    theVFName = ['VF' cfdGetFluidTagUsingIndex(iFluid)];
    VF=cfdGetMeshField(theVFName);
    sumVF=sumVF+VF.phi;
end
