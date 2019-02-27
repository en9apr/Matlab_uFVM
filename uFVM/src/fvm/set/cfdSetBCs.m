function theEquation = cfdSetBCs(theName,theBCCells)
%
%
%

theEquationName = cfdConvertName(theName);

theEquation = cfdGetField(theEquationName);


theMesh = cfdGetMesh;
theEquation.boundaryConditions  = cfdSetupBoundaryConditions(theMesh,theBCCells);



cfdSetField(theEquation);
