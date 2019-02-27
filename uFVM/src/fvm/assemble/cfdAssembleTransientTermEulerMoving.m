function theFluxes = cfdAssembleTransientTermEulerMoving(theEquationName,theTerm,iComponent)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================
%---------------------- Start Term Info ---------------------
%
theMesh = cfdGetMesh;

theEquation = cfdGetModel(theEquationName);
theFluidTag = cfdGetFluidTag(theEquationName);
theRhoName = theEquation.rhoName;
dt = cfdGetDt;

iElements = 1:theMesh.numberOfElements;

theFluxes.FLUXCE(iElements)    = 0;
theFluxes.FLUXCEOLD(iElements) = 0;
theFluxes.FLUXTE(iElements)    = 0;
%
%---------------------- End Term Info ----------------------
%
theEquationMeshField = cfdGetMeshField(theEquationName);
phi = theEquationMeshField.phi(iElements,iComponent);
theEquationMeshFieldOld = cfdGetMeshField(theEquationName,'Elements','Step1');
phi_old = theEquationMeshFieldOld.phi(iElements,iComponent);
%
% Temperature Equation
%
%  
theEquationBaseName = cfdGetBaseName(theEquationName);

if(strcmp(theEquationBaseName,'Temperature'))
    %
    theDensityField = cfdGetMeshField(['Density' theFluidTag]);    
    theSpecificHeatField = cfdGetMeshField(['SpecificHeat' theFluidTag]);
    density = theDensityField.phi(iElements);    
    cp = theSpecificHeatField.phi(iElements);
    
    
    theDensityFieldOld = cfdGetMeshField(['Density' theFluidTag],'Elements','Step1');    
    theSpecificHeatFieldOld = cfdGetMeshField(['SpecificHeat' theFluidTag],'Elements','Step1');
    density_old = theDensityFieldOld.phi(iElements);
    cp_old= theSpecificHeatFieldOld.phi(iElements);

    volumes = [theMesh.elements(iElements).volume]';
    
    theFluxes.FLUXCE(iElements)    =   volumes .* cp .* density / dt;
    theFluxes.FLUXCEOLD(iElements) = - volumes .* cp_old .* density_old / dt;
    theFluxes.FLUXTE(iElements)    =  theFluxes.FLUXCE .* phi + theFluxes.FLUXCEOLD .* phi_old ;
%
% VF Equations
%
elseif(strcmp(theEquationBaseName,'VF')==true)
    %
    theDensityField = cfdGetMeshField(['Density' theFluidTag]);    
    density = theDensityField.phi(iElements);
    density_old = theDensityField.phi_old(iElements);
    size(density)
    size(volumes)
    
    volumes = [theMesh.elements(iElements).volume]';
    
    theFluxes.FLUXCE(iElements)    =   volumes .* density / dt;
    theFluxes.FLUXCEOLD(iElements) = - volumes .* density_old / dt;
    theFluxes.FLUXTE(iElements)    = theFluxes.FLUXCE .* phi + theFluxes.FLUXCEOLD .* phi_old ;

%
% All Other Equations 
%
else
    %
    % specifiy the Term Coefficient Field
    %
    theRhoMeshField = cfdGetMeshField(theRhoName); 
    
    if(isempty(theRhoMeshField))
        theFluxes.FLUXCE(iElements) = 0;
        theFluxes.FLUXCEOLD(iElements) = 0;
        theFluxes.FLUXVE(iElements) = 0;
        theFluxes.FLUXTE(iElements) = 0;
    else 
        rho = theRhoMeshField.phi(iElements);

        theRhoMeshFieldOld = cfdGetMeshField(theRhoName,'Elements','Step1');    
        rho_old = theRhoMeshFieldOld.phi(iElements);


        volumes = [theMesh.elements(iElements).volume]';
        oldvolumes= [theMesh.elements(iElements).OldVolume]';
 
        theFluxes.FLUXCE1(iElements) =  volumes .* rho / dt ;
        theFluxes.FLUXCE2(iElements) = -oldvolumes .* rho / dt;        
        theFluxes.FLUXCE(iElements) = theFluxes.FLUXCE1(iElements)+theFluxes.FLUXCE1(iElements)+ theFluxes.FLUXCE2(iElements) ;       
        theFluxes.FLUXCEOLD(iElements) = - volumes .* rho_old / dt;
        
        theFluxes.FLUXTE(iElements)    = theFluxes.FLUXCE .* phi' + theFluxes.FLUXCEOLD  .* phi_old';
    end
end

    
