function fieldExists = cfdCheckIfFieldExists(theFieldName)
%===================================================

%  written by the CFD Group @ AUB, Fall 2017
%===================================================
global Domain;

fieldExists = false;
theNumberOfFields = length(Domain.fields);
for iField=1:theNumberOfFields
    theName = Domain.fields{iField}.name;
    if(strcmp(theName,theFieldName))
        fieldExists = true;
    end
end

