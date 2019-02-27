function cfdSetModel(theField)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================

global Domain;

theFieldIndex = 0;
theNumberOfFields = length(Domain.fields);
for iField=1:theNumberOfFields
    if(strcmp(theField.name,Domain.fields{iField}.name))
        Domain.fields{iField} = theField;
        theFieldIndex = iField;
    end
end
if(theFieldIndex==0)
    Domain.fields{theNumberOfFields+1} = theField;
end
end