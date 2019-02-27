function cfdSetupProperties
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2017
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function set up properties
%--------------------------------------------------------------------------

% Retrieve transport properties from data base
global Domain;
transportProperties = Domain.foam.transportProperties;

% Check Class of Application
applicationClass = cfdGetApplicationClass;
if(strcmp(applicationClass, 'incompressible') || strcmp(applicationClass, 'basic'))
    propertyNames = fieldnames(transportProperties);
    for iPropertyName=1:length(propertyNames)
        propertyName = propertyNames{iPropertyName};
        property = getfield(transportProperties, propertyName);
        propertyValue = num2str(property{2,2});
        cfdSetupProperty(propertyName, 'constant', propertyValue);
    end
elseif(strcmp(applicationClass, 'compressible'))
    
elseif(strcmp(applicationClass, 'multiphase'))
    
elseif(strcmp(applicationClass, 'heatTransfer'))
    
end





