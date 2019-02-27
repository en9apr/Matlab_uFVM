function cfdUpdateTransientProperty(thePropertyName)

theMeshField = cfdGetMeshField(thePropertyName,'Elements','Step0');
cfdSetMeshField(theMeshField,'Step1')