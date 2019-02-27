function cfdUpdateTransientTermCoefficients(theEquationName)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function stores current time step fields to previous ones
%--------------------------------------------------------------------------

% Update terms coefficients (excluding p equation)
if strcmp(theEquationName, 'U')
    theEquation = cfdGetModel(theEquationName);
    theTerms = theEquation.terms;
    for iTerm=1:length(theTerms)
        theTerm = theTerms{iTerm};
        if strcmp(theTerm.fvmOperator, 'div') && ~strcmp(theEquationName, 'U')
            theMeshField = cfdGetMeshField(['psi_',theEquationName,'eq'],'Faces','Step0');
            cfdSetMeshField(theMeshField,'Step1');
        elseif strcmp(theTerm.fvmOperator, 'laplacian')
            theMeshField = cfdGetMeshField(['gamma_',theEquationName,'eq'],'Faces','Step0');
            cfdSetMeshField(theMeshField,'Step1');
        end
    end
    
    % Always update the rho field
    theMeshField = cfdGetMeshField(['rho_',theEquationName,'eq'],'Elements','Step0');
    cfdSetMeshField(theMeshField,'Step1');
elseif ~strcmp(theEquationName, 'p')
    theEquation = cfdGetModel(theEquationName);
    theTerms = theEquation.terms;
    for iTerm=1:length(theTerms)
        theTerm = theTerms{iTerm};
        if strcmp(theTerm.fvmOperator, 'div') && ~strcmp(theEquationName, 'U')
            theMeshField = cfdGetMeshField(['psi_',theEquationName,'eq'],'Faces','Step0');
            cfdSetMeshField(theMeshField,'Step1');
        elseif strcmp(theTerm.fvmOperator, 'laplacian')
            theMeshField = cfdGetMeshField(['gamma_',theEquationName,'eq'],'Faces','Step0');
            cfdSetMeshField(theMeshField,'Step1');
        end
    end
    
    % Always update the rho field
    theMeshField = cfdGetMeshField(['rho_',theEquationName,'eq'],'Elements','Step0');
    if ~isempty(theMeshField)
        cfdSetMeshField(theMeshField,'Step1');   
    end
end

%
% update mdot_f
%
if strcmp(theEquationName, 'U')
    applicationClass = cfdGetApplicationClass;
    if strcmp(applicationClass, 'multiphase')
        theNumberOfFluids = cfdGetNumberOfFluids;
        for iFluid=1:theNumberOfFluids
            cfdUpdateTransientMdot(['mdot_f', num2str(iFluid)]);
        end
    end
    cfdUpdateTransientMdot('mdot_f');
end
