function theField = cfdGetModel(theFieldName)


theField={};
theFields = cfdGetModels;
theNumberOfFields = cfdGetNumberOfFields;
for iField=1:theNumberOfFields
    if(strcmp(theFieldName,theFields{iField}.name) || strcmp(theFieldName,theFields{iField}.userName))
        theField = theFields{iField};
        return
    end
end
