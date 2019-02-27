function cfdInitializeConstantFields
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function initializes the constant fields
%--------------------------------------------------------------------------

global Domain;

fieldNames = fieldnames(Domain.foam.fields);
theNumberOfFieldNames = length(fieldNames);

theEquationNames = cfdGetEquationNames;
theNumberOfEquations = length(theEquationNames);

for iField=1:theNumberOfFieldNames
    theFieldName = fieldNames{iField};    
    
    % Skip mdot_f field if found. mdot_f is always available in the time
    % directories except for 0 directory.
    if strcmp(theFieldName, 'mdot_f')
        continue;
    end
    
    isConstant = true;
    for iEquation=1:theNumberOfEquations
        theEquationName = theEquationNames{iEquation};
        if strcmp(theEquationName, theFieldName)
            isConstant = false;
            break;
        end
    end
    
    if isConstant
        theField = Domain.foam.fields.(theFieldName);
        
        % Create the mesh field corresponding to the constant field
        if strcmp(theField.class, 'volScalarField')
            cfdSetupMeshField(theFieldName,'Elements','Scalar');
        else
            cfdSetupMeshField(theFieldName,'Elements','Vector');            
        end
        cfdCreateMeshField(theFieldName, theField);
    end
end