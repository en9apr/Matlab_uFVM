function cfdCorrectEquation(theEquationName,iComponent)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function corrects the equations
%--------------------------------------------------------------------------

theEquationBaseName = cfdGetBaseName(theEquationName);
theEquation = cfdGetModel(theEquationName);
theEquationType = theEquation.type;
%
% VELOCITY Correction
%
if strcmp(theEquationType,'Vector')
    if(strcmp(theEquationBaseName,'U'))
        cfdCorrectVelocityEquation(theEquationName,iComponent);
    else
        cfdCorrectScalarEquation(theEquationName,iComponent);
    end
    %
    % PRESSURE Correction
    %
elseif strcmp(theEquationBaseName,'p')
    cfdCorrectPPField;
    
    % Correct U
    cfdCorrectVelocityField;
    
    % Correct p
    cfdCorrectPressureEquation;
    
    % Correct rho, if the flow is compressible
    applicationClass = cfdGetApplicationClass;
    if strcmp(applicationClass, 'compressible')
        cfdCorrectDensityField;
    end
    
    % Correct mdot_f
    applicationClass = cfdGetApplicationClass;
    if strcmp(applicationClass, 'multiphase')
        cfdCorrectMixtureMdotField;
    else
        cfdCorrectMdotField;
    end
       
else
    cfdCorrectScalarEquation(theEquationName,iComponent);
end


end
