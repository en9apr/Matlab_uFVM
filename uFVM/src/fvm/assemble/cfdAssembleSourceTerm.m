function cfdAssembleSourceTerm(theTerm,iComponent)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function assembles the source term
%--------------------------------------------------------------------------

% Frequent source terms in the momentum and energy equation. Treating each
% of these in a special way is better than direct evaluating as it consumes
% time
if strcmp(theTerm.formula, 'grad(p)')
    cfdAssemblePressureGradientTerm(theTerm,iComponent);
elseif strcmp(theTerm.formula, 'div(mu*transp(grad(U)))') || strcmp(theTerm.formula, 'div(nu*transp(grad(U)))')
    cfdAssembleExplicitStressTerm(theTerm,iComponent); 
elseif strcmp(theTerm.formula, '2/3*grad(mu*div(U))')
    cfdAssembleBulkViscosityTerm(theTerm,iComponent);
elseif strcmp(theTerm.formula, 'rho*g') || strcmp(theTerm.formula, 'g')
    cfdAssembleGravitationalForceTerm(theTerm,iComponent);    
elseif strcmp(theTerm.formula, 'rho*g*beta*(T-TRef)') || strcmp(theTerm.formula, 'g*beta*(T-TRef)')
    cfdAssembleBuoyancyTerm(theTerm,iComponent);
elseif strcmp(theTerm.formula, 'rho*T*DDt(Cp)')
    cfdAssembleSpecificHeatTerm(theTerm);
elseif strcmp(theTerm.formula, 'DDt(p)')
    cfdAssembleSubstantialDerivativeTerm(theTerm);
elseif strcmp(theTerm.formula, '2/3*mu*PSI')
    cfdAssembleDissipationTerm(theTerm);
elseif strcmp(theTerm.formula, 'mu*PHI')
    cfdAssembleViscousDissipationTerm(theTerm);
else
    % If the source term is not recognized by the above terms, it is
    % calculated by evaluating the terms directly
    
    % Get mesh info
    theMesh = cfdGetMesh;
    iElements = 1:theMesh.numberOfElements;
    volume = [theMesh.elements(iElements).volume]';
    
    % Initialize volume fluxes
    theFluxes.FLUXCE(iElements, 1) = 0;
    theFluxes.FLUXCEOLD(iElements, 1) = 0;
    theFluxes.FLUXVE(iElements, 1) = 0;
    theFluxes.FLUXTE(iElements, 1) = 0;
    theFluxes.FLUXTEOLD(iElements, 1) = 0;
    
    % Identify source term. Some source terms are unique or frequent like field
    % gradients (i.e. grad(p), div(mu*transp(grad(U))), 2/3*grad(mu*div(U))).
    % If so, then it is easier since field gradient are available in hand and
    % no need to re-compute it. Other source terms may be treated in a quicker
    % way
    theFormula = theTerm.formula;       
    
    % If the formula is grad(phi), where phi is any of the equation field names
    theEquationNames = cfdGetEquationNames;
    for iEquation=1:length(theEquationNames)
        if strcmp(theFormula, ['grad(', theEquationNames{iEquation},')'])
            % Get the stored gradient of the mesh field
            theMeshField = cfdGetMeshField(theEquationNames{iEquation});
            S = theMeshField.phiGradient(iElements, iComponent);
            
            % Assemble Source Term as element flux
            theFluxes.FLUXTE(iElements) = theTerm.sign * S .* volume;
            cfdAssembleIntoGlobalMatrixElementFluxes(theFluxes);
            return;
        end
    end    
    
    % If source term is not standard
    S = cfdEvaluateNonstandardSourceTerm(theTerm, iComponent);
    
    % Assemble Source Term as element flux
    theFluxes.FLUXTE(iElements) = theTerm.sign * S .* volume;
    cfdAssembleIntoGlobalMatrixElementFluxes(theFluxes);
end

end