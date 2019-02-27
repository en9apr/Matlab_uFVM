function cfdUpdateDensityField(theDensityModel)

theEquation = cfdGetModel(theDensityModel);
if(isfield(theEquation,'model'))
    theModelName = theEquation.model;
    theFluidTag = theEquation.tag;
    
    thePressureField = cfdGetMeshField('Pressure');
    P = thePressureField.phi;
    
    theTemperatureField = cfdGetMeshField(['Temperature' theFluidTag]);
    T = theTemperatureField.phi;
    
    theGasConstantField = cfdGetMeshField(['GasConstant' theFluidTag]);
    R = theGasConstantField.phi;
    
    theDensityField = cfdGetMeshField(['Density' theFluidTag]);
    
    if strcmp(theModelName,'Ideal Gas')
        rho = P./(R.*T);
        theDensityField.phi = rho;
        cfdSetMeshField(theDensityField);
    elseif strcmp(theModelName,'Real Gas')
        
    end
end
