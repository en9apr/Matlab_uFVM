function cfdSetGradient(theName,theGradientName)



theEquationName = cfdConvertName(theName);

theScalarField = cfdGetModel(theEquationName);

theScalarField.gradientName = theGradientName;

cfdSetMeshField(theScalarField);