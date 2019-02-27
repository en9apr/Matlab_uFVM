function [rmsResidual, maxResidual] = cfdAssembleEquationTerms(theEquationName,iComponent)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017 
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function assembles equation terms
%--------------------------------------------------------------------------

if(nargin==1)
    iComponent = 1;
end

% get theScalarName
theEquation = cfdGetModel(theEquationName);

% check if equation is to be assembled
theNumberOfTerms = length(theEquation.terms);

% Assemble Coefficients
for iTerm = 1:theNumberOfTerms
    theTerm = theEquation.terms{iTerm};
    if strcmp(theTerm.name,'Transient')
        cfdAssembleTransientTerm(theEquationName,theTerm,iComponent);
    elseif strcmp(theTerm.name,'Convection')
        cfdAssembleConvectionTerm(theEquationName,theTerm,iComponent);
    elseif strcmp(theTerm.name,'Diffusion')
        cfdAssembleDiffusionTerm(theEquationName,theTerm,iComponent);
    elseif strcmp(theTerm.name,'Source')
        cfdAssembleSourceTerm(theTerm,iComponent); 
    elseif strcmp(theTerm.name,'Implicit Source')        
        cfdAssembleSourceTerm(theTerm,iComponent); % not implemented yet
    elseif strcmp(theTerm.name,'mdot_f')
        cfdAssembleMdotTerm(theEquationName,theTerm);  
    elseif strcmp(theTerm.name,'mdot_mixture_f')
        cfdAssembleMixtureMdotTerm(theEquationName,theTerm);          
    else
        error('\n%s\n',[theTerm.name,' term is not defined']);
    end
end

% Add false transience to equations other than pressure equation if and
% only if the simulation is steady
if ~cfdIsTransient && ~strcmp(theEquationName,'p')
    cfdAssembleFalseTransientTerm(theEquationName,iComponent);
end

% Apply under-relaxation for all equations other than pressure equation
if ~strcmp(theEquationName,'p')
    cfdApplyURF(theEquationName);
end

% Compute RMS and MAX Residuals
[rmsResidual, maxResidual] = cfdComputeNormalizedResidual(theEquationName, iComponent);

end
