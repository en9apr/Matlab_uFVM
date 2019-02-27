function cfdSetPhysicalBCForPatch(iPatch,physicalBC)

global Domain;

Domain.mesh.boundaries(iPatch).type=physicalBC;