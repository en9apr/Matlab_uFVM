function cfdUpdateTransientEquation(theEquationName)

theMeshField = cfdGetMeshField(theEquationName,'Elements','Step0');
cfdSetMeshField(theMeshField,'Step1')