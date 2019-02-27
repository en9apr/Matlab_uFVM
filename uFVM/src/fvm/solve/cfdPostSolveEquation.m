function cfdPostSolveEquation(theEquationName,iComponent)
%===================================================
%
%  written by the CFD Group @ AUB, Fall 2006
%===================================================
if(nargin==1)
    iComponent=1;
end

theEquationBaseName=cfdGetBaseName(theEquationName);
theFluidTag=cfdGetFluidTag(theEquationName);
theFluidIndex = cfdGetFluidIndex(theEquationName);
%
%
% Temperature Equations
%
if(strcmp(theEquationBaseName,'T'))
    %
    % fvmUpdateVolumeFractionField
    %
end

if(strcmp(theEquationBaseName,'H'))
    
    %
    % Updating Temperature
    
    thePropertyFieldName = ['Temperature' theFluidTag];
    thePropertyField = cfdComputePropertyField(thePropertyFieldName);
    
    if (isfield(thePropertyField, {['phiGradient' theFluidTag]}))
        thePropertyField = cfdUpdateGradient(thePropertyFieldName);
    end
    
    cfdUpdateField(thePropertyField);
    %
    % Updating Liquid Concentration fraction
    thePropertyField = cfdGetModel(['Cl' theFluidTag]);
    thePropertyField = cfdComputePropertyField(thePropertyField);
    
    if (isfield(thePropertyField, {['phiGradient' theFluidTag]}))
        thePropertyField = fvmUpdateGradient(thePropertyField);
    end
    
    cfdUpdateField(thePropertyField);
    
end


