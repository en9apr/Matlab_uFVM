function cfdWriteResults(currentGlobalIteration, currentTime)
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function writes the results at each write interval
%--------------------------------------------------------------------------

global Domain;

% Get control settings
controlDict = Domain.foam.controlDict;
writeControl = controlDict{7, 2};
writeInterval = controlDict{8, 2};

% Get equation names
theEquationNames = cfdGetEquationNames;

if strcmp(writeControl, 'timeStep')
    if rem(currentGlobalIteration, writeInterval)==0
        mkdir(num2str(currentTime));
        for iEquation=1:length(theEquationNames)
            theEquationName = theEquationNames{iEquation};
            theField = cfdGetMeshField(theEquationName);
            cfdWriteOpenFoamField(theField, currentTime);
            
            % Write mdot_f field if U equation is being solved
            if strcmp(theEquationName, 'U')
                theField = cfdGetMeshField('mdot_f', 'Faces');
                cfdWriteOpenFoamField(theField, currentTime);
            end
        end
    end
end
