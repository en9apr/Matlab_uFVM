function cfdSetupSpecies(theSpecieUserName)


theFluid = cfdGetFluidUsingName(theSpecieUserName);
theFluidTag = theFluid.tag;
theFluidUserName = theFluid.userName;

theNumberOfSpecies = length(theFluid.species);
theSpecieIndex = theNumberOfSpecies+1;

theSpecie.userName = theSpecieUserName;
theSpecie.name = theSpecieUserName;
theSpecie.index = theSpecieIndex;
theSpecie.tag = ['_specie0' num2str(theSpecieIndex) theFluidTag];
theSpecie.type = ''; 

theFluid{theSpecie.index}=theSpecie;
cfdSetFluid(theFluid);
