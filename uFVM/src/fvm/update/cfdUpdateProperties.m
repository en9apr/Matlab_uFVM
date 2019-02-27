function cfdUpdateProperties
%===================================================

%  written by the CFD Group @ AUB, Fall 2017
%===================================================
%
% This function updates the properties defined by the user. The property
% could be a constant, a formula or a model.
%
thePropertyNames = cfdGetPropertyNames;
theNumberOfProperties = length(thePropertyNames);

for iProperty = 1:theNumberOfProperties
    thePropertyName = thePropertyNames{iProperty};
    cfdUpdateProperty(thePropertyName);
    
    theProperty = cfdGetModel(thePropertyName);
    % If the property includes a gradient to be computed, then it is done
    % here
    if (isfield(theProperty, {'phiGradient'}))
        cfdUpdateGradient(thePropertyName); % this part has been added
    end
    
end
