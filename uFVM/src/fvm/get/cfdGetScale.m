function phiScale = cfdGetScale(theEquationUserName)

theField = cfdGetMeshField(theEquationUserName);
phiScale = theField.phiScale;