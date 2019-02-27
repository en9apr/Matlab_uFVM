function cfdUpdateTransientMdotFields(theEquationName)

theMeshField = cfdGetMeshField(theEquationName,'Step0');
cfdSetMeshField(theMeshField,'Step1')