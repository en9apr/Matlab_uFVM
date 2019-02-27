function theSourceFieldNames = cfdGetSourceFieldNames;
%===================================================

%  written by the CFD Group @ AUB, Fall 2006
%===================================================

global Domain;

theSourceFieldNames = {};
theNumberOfFields = length(Domain.fields);
theNumberOfSourceFields = 0;
for iField=1:theNumberOfFields
    if(strcmp(Domain.fields{iField}.type,'Source'))
       theNumberOfSourceFields = theNumberOfSourceFields+1;
       theSourceFieldNames{theNumberOfSourceFields} = Domain.fields{iField}.name; 
    end
end