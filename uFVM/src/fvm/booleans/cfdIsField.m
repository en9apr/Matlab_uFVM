function isField = cfdIsField(theFullFieldName)
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================

global Domain;

isField = false;

theFieldName = cfdConvertFormula(theFullFieldName);
% theField = 0;
theField = {};
theNumberOfFields = length(Domain.fields);
for iField=1:theNumberOfFields
    if(strcmp(theFieldName,Domain.fields{iField}.name))
        isField = true;
        return
    end
end
