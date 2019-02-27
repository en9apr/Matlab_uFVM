function cfdSetSpecies(theSpecie)


theFluid{theSpecie.index} = theSpecie;
cfdUpdateFluid(theFluid);
