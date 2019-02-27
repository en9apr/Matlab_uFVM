function cfdGetUpwindElementIndices(theFluidTag)


if(nargin<1)
    mdotName = ['Mdot' theFluidTag];
end

mdot = cfdGetMeshField(mdotName);

