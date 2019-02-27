function bodyForces = cfdGetBodyForces
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017 
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function collects the terms from the momentum equation that are 
%   considered body forces like gravitational force or electrostatic force
%--------------------------------------------------------------------------

bodyForces = {};
theVelocityEquation = cfdGetModel('U');
theTerms = theVelocityEquation.terms;
iBodyForce = 1;
for iTerm=1:length(theTerms)
    theTerm = theTerms{iTerm};
    if strcmp(theTerm.name, 'Source')
        % Check if the term includes a gravity field
        fieldNames = indentifyFields(theTerm.formula);
        if any(ismember(fieldNames, 'g'))
            bodyForces{iBodyForce, 1} = theTerm.formula;
            iBodyForce = iBodyForce + 1;
        end        
    end
end

