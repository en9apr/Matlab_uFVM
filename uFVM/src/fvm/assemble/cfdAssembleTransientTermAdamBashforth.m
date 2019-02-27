function theCoefficients = cfdAssembleTransientTermAdamBashforth(theCoefficients,theScalarFieldName,dt,theTerm)
%===================================================
%
% This function computed the correction to the UPWIND convection assembly
% to yield the SOU scheme assemble using the deferred correction method
%
%  written by the CFD Group @ AUB, Fall 2006
%===================================================

theFluidIndex = cfdGetFluidIndex(theScalarFieldName);
theFluidTag = cfdGetFluidTag(theScalarFieldName);
theFieldBaseName = cfdGetBaseName(theScalarFieldName);
theScalarField = cfdGetField(theScalarFieldName);
% Account for Term Type
%
theTermType = theTerm.type;
if(strcmp(theTermType,'Residual'))
  cLHS=1;
  cRHS=1;
elseif(strcmp(theTermType,'False'))
    cLHS=1;
    cRHS=0;
elseif(strcmp(theTermType,'Deferred'))
    cLHS=0;
    cRHS=1;
else
    error([theTerm.name '<<<< has not type defined'])
end
%
% specifiy the Term Field
%
if(isempty(theTerm.variable))
   theTermField = theScalarField;
else
   theTermField = cfdGetField(theTerm.variable);
end
%
% specifiy the Term Coefficient Field
%
if(isempty(theTerm.coefficient))
   theTermCoefficientField = cfdGetField(theScalarField.rho);
else
   theTermCoefficientField = cfdGetField(theTerm.coefficient);
end
%---------------------- End Term Info ----------------------
a = theCoefficients.LHS;
b = theCoefficients.RHS;

theMesh = cfdGetMesh;
theElements = theMesh.elements;
theFaces = theMesh.faces;


rho = theTermCoefficientField.phi;
rho_old = theTermCoefficientField.phiOld;
rho_oldold = theTermCoefficientField.phiOldOld;
%
phi = theScalarField.phi;
phiOld = theScalarField.phiOld;
phiOldOld = theScalarField.phiOldOld;


theVFField = cfdGetField(['VF' theFluidTag]);
vf = theVFField.phi;
%   

%//////////////////////////////////////////////////////////
%
% Loop and assemble the corrections
%

theNumberOfElements = theMesh.numberOfElements;
if(strcmp(theFieldBaseName,'Temperature'))
    %
    theDensityField = cfdGetField(['Density' theFluidTag]);    
    theSpecificHeatField = cfdGetField(['SpecificHeat' theFluidTag]);

    density = theDensityField.phi;
    density_old = theDensityField.phiOld;
    density_oldold = theDensityField.phiOldOld;

    
    cp = theSpecificHeatField.phi;
    cp_old= theSpecificHeatField.phiOld;
    cp_oldold= theSpecificHeatField.phiOldOld;
    
    for iElement = 1:theNumberOfElements
        volume = theElements(iElement).volume;
      % Compute Fluxes     
        %
        FLUXCE =     1.5*volume*cp(iElement)*density(iElement)/dt;
        FLUXCE_OLD = -2*volume*cp_old(iElement)*density_old(iElement)/dt;
        FLUXCE_OLDOLD= 0.5*volume*cp_oldold(iElement)*density_oldold(iElement)/dt;
        FLUXVE = 0;
        FLUXTE =   FLUXCE*phi(iElement) + FLUXCE_OLD*phiOld(iElement) + FLUXCE_OLDOLD*phiOldOld(iElement)+  FLUXVE;
        %
        % Assemble Term
        %
        a(iElement,iElement) = a(iElement,iElement) + cLHS*FLUXCE*vf(iElement);
        b(iElement) = b(iElement) - cRHS*FLUXTE*vf(iElement) ;
    end
  elseif(strcmp(theFieldBaseName,'VF')==true)
  %
  
  %
  %
  
    for iElement = 1:theNumberOfElements
       volume = theElements(iElement).volume;
        %
        % Compute Fluxes     
        %
        FLUXCE =     1.5*volume*rho(iElement)/dt;
        FLUXCE_OLD = -2*volume*rho_old(iElement)/dt;
        FLUXCE_OLDOLD = 0.5*volume*rho_oldold(iElement)/dt;
        FLUXVE = 0;
        FLUXTE =   FLUXCE*phi(iElement) + FLUXCE_OLD*phiOld(iElement) + FLUXCE_OLDOLD*phiOldOld(iElement)+  FLUXVE;
        %
        % Assemble Term
        %
        a(iElement,iElement) = a(iElement,iElement) + cLHS*FLUXCE;
        b(iElement) = b(iElement) - cRHS*FLUXTE;
    
    end
  else
    for iElement = 1:theNumberOfElements
        volume = theElements(iElement).volume;
        %
        % Compute Fluxes     
        %
        FLUXCE =      1.5*volume*rho(iElement)/dt;
        FLUXCE_OLD = -2*volume*rho_old(iElement)/dt;
        FLUXCE_OLDOLD = 0.5*volume*rho_oldold(iElement)/dt;
        FLUXVE = 0;
        FLUXTE =   FLUXCE*phi(iElement) + FLUXCE_OLD * phiOld(iElement)+ FLUXCE_OLDOLD *phiOldOld(iElement) + FLUXVE;
        %
        % Assemble Term
        %
        a(iElement,iElement) = a(iElement,iElement) + cLHS*FLUXCE*vf(iElement);
        b(iElement) = b(iElement) - cRHS*FLUXTE*vf(iElement);
        
        
    end
end
  
    
   % b(iElement) = b(iElement) + a(iElement,iElement);
   % a(iElement,iElement) = a(iElement,iElement) + FLUXTE;
   % b(iElement2) = b(iElement2) - FLUXTf*vf_f(iFace);

theCoefficients.LHS = a;
theCoefficients.RHS = b;
